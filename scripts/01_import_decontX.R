rm(list=ls());gc()
setwd("/home/carinotus1/课题/202603新生期口服过敏原免疫细胞图谱")
getwd()

#!/usr/bin/env Rscript
# Step 01: import + QC (per batch or per sample)
# Input: metadata/samples.csv + 10x matrices
# Output: QC plots + milestone Seurat objects in data/processed/
message("[01] import + decontX: start")


#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(celda) # decontX
  library(SingleCellExperiment)
})

meta <- read.csv("metadata/samples_auto.csv", stringsAsFactors = FALSE);view(meta)

dir.create("results/figures/01_qc_diag", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/01_preQC", recursive = TRUE, showWarnings = FALSE)
dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)

# 粗阈值
prefilter <- list(min_feat = 100, min_count = 200, max_mt = 40)
decont <- list(res = 0.6, pcs = 1:20)




# 确认基因格式，小鼠都是mt-和Rp[sl]
detect_patterns <- function(genes){
  mt_pat <- if (any(grepl("^mt-", genes))) "^mt-" else if (any(grepl("^MT-", genes))) "^MT-" else "^mt-"
  ribo_pat <- if (any(grepl("^Rp[sl]", genes))) "^Rp[sl]" else if (any(grepl("^RP[SL]", genes))) "^RP[SL]" else "^Rp[sl]"
  list(mt=mt_pat, ribo=ribo_pat)
}


# function process_one：批量读入+deContX去除环境污染
process_one <- function(sample_id, mrna_path){
  
  message("[", sample_id, "] read")
  counts <- Read10X(data.dir = mrna_path)
  
  obj <- CreateSeuratObject(counts, project = sample_id, min.cells = 3, min.features = 200)
  obj$orig.ident <- sample_id
  
  # mt/ribo pattern（mouse优先，小写/大写自动兼容）
  genes <- rownames(obj)
  mt_pat   <- if (any(grepl("^mt-", genes))) "^mt-" else if (any(grepl("^MT-", genes))) "^MT-" else "^mt-"
  ribo_pat <- if (any(grepl("^Rp[sl]", genes))) "^Rp[sl]" else if (any(grepl("^RP[SL]", genes))) "^RP[SL]" else "^Rp[sl"
  
  obj$mt_percent   <- PercentageFeatureSet(obj, pattern = mt_pat)
  obj$ribo_percent <- PercentageFeatureSet(obj, pattern = ribo_pat)
  
  p1 <- VlnPlot(
    obj, layer = "counts",
    features = c("nCount_RNA","nFeature_RNA","mt_percent","ribo_percent"),
    pt.size = 0, ncol = 4
  ) + plot_annotation(title = paste0(sample_id, " raw"))
  
  # 轻预过滤（给decontX用）
  obj0 <- subset(obj, subset = nFeature_RNA > 100 & nCount_RNA > 200 & mt_percent < 40)
  
  message("[", sample_id, "] rough cluster for decontX")
  obj0 <- NormalizeData(obj0, verbose = FALSE)
  obj0 <- FindVariableFeatures(obj0, nfeatures = 2000, verbose = FALSE)
  obj0 <- ScaleData(obj0, verbose = FALSE)
  obj0 <- RunPCA(obj0, npcs = 30, verbose = FALSE)
  obj0 <- FindNeighbors(obj0, dims = 1:20, verbose = FALSE)
  obj0 <- FindClusters(obj0, resolution = 0.6, verbose = FALSE)
  
  message("[", sample_id, "] decontX")
  sce <- as.SingleCellExperiment(obj0)
  sce <- celda::decontX(sce, z = obj0$seurat_clusters)
  
  contam <- colData(sce)$decontX_contamination
  dc_counts <- decontXcounts(sce)
  
  # 兜底
  contam[is.na(contam)] <- 0
  
  message("[", sample_id, "] rebuild seurat (decontX counts)")
  obj_dc <- CreateSeuratObject(dc_counts, project = sample_id, min.cells = 3, min.features = 200)
  obj_dc$orig.ident <- sample_id
  obj_dc$decontX_contam <- contam
  
  genes2 <- rownames(obj_dc)
  mt_pat2   <- if (any(grepl("^mt-", genes2))) "^mt-" else if (any(grepl("^MT-", genes2))) "^MT-" else "^mt-"
  ribo_pat2 <- if (any(grepl("^Rp[sl]", genes2))) "^Rp[sl]" else if (any(grepl("^RP[SL]", genes2))) "^RP[SL]" else "^Rp[sl]"
  obj_dc$mt_percent   <- PercentageFeatureSet(obj_dc, pattern = mt_pat2)
  obj_dc$ribo_percent <- PercentageFeatureSet(obj_dc, pattern = ribo_pat2)
  
  p2 <- VlnPlot(
    obj_dc, layer = "counts",
    features = c("nCount_RNA","nFeature_RNA","mt_percent","ribo_percent","decontX_contam"),
    pt.size = 0, ncol = 5
  ) + plot_annotation(title = paste0(sample_id, " decontX"))
  
  p3a <- FeatureScatter(obj_dc, "nCount_RNA", "nFeature_RNA")
  p3b <- FeatureScatter(obj_dc, "nCount_RNA", "mt_percent")
  p3  <- (p3a + p3b) + plot_layout(ncol = 2) + plot_annotation(title = paste0(sample_id, " scatter (pre strict QC)"))
  
  pdf(file = file.path("results/figures/01_qc_diag", paste0(sample_id, "_diag.pdf")), width = 12, height = 15)
  print(p1 / p2 / p3)
  dev.off()
  
  p4 = p1 / p2 / p3
  ggsave(file = file.path("results/figures/01_qc_diag", paste0(sample_id, "_diag.png")),plot = p4, width = 12, height = 15)
  
  saveRDS(obj_dc, file = file.path("data/processed/01_preQC", paste0(sample_id, ".rds")))
  invisible(obj_dc)
}

# test one sample
# process_one("CONT_D_4", meta$mrna_path[meta$sample_id=="CONT_D_4"])
# [CONT_D_4] rebuild seurat (decontX counts)
# done

# run in bash
for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  mp  <- meta$mrna_path[i]
  message("[01] ", sid)
  
  tryCatch({
    process_one(sid, mp)
  }, error = function(e){
    msg <- paste0(Sys.time(), "  ", sid, "  ERROR: ", e$message)
    write(msg, file = "results/logs/01_import_decontX_errors.txt", append = TRUE)
    message(msg)
  })
}
