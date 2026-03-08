#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(DoubletFinder)
  library(readr)
})

# -----------------------------
# user settings
# -----------------------------
in_dir  <- "data/processed/02_QC"
out_dir <- "data/processed/03_doublet"
fig_dir <- "results/figures/03_doublet"
tab_dir <- "results/tables"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# 只标注 doublet，先不删
remove_doublets <- FALSE   # TRUE：输出去除后对象；FALSE：只标注

# PCs / 聚类参数（DoubletFinder也用这套）
npcs <- 50
PCs_use <- 1:15
cluster_res_for_homotypic <- 0.3

# pK 策略：
# - "per_sample": 每样本都 paramSweep 找 pK（慢但精确）
# - "reuse_one":  只对一个代表样本找 pK，其他样本复用（快，推荐）
pk_mode <- "reuse_one"

# 若 reuse_one：用哪个样本来找 pK？
# （建议选细胞数中等、质量不错的）
pk_reference_sample <- NULL  # 例如 "CONT_C_1"；NULL 则自动选中位数细胞量的样本

# Doublet rate 经验估计（按细胞数粗估）
estimate_doublet_rate <- function(nCells){
  dplyr::case_when(
    nCells < 5000  ~ 0.04,
    nCells < 10000 ~ 0.06,
    nCells < 20000 ~ 0.075,
    nCells < 40000 ~ 0.10,
    TRUE           ~ 0.12
  )
}

# -----------------------------
# helpers
# -----------------------------
stop_msg <- function(...) stop(paste0(...), call. = FALSE)

prep_for_df <- function(obj){
  # normalize (LogNormalize) — 你之前就是这条路线
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = PCs_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_res_for_homotypic, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = PCs_use, verbose = FALSE)
  obj
}

find_best_pk <- function(obj, sample_id){
  sweep.res.list <- paramSweep(obj, PCs = PCs_use, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # 保存 BCmetric 曲线
  pdf(file.path(fig_dir, paste0(sample_id, "_findpK.pdf")), width = 7, height = 4)
  print(ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_line() + geom_point())
  dev.off()
  
  mpK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  mpK
}

run_doubletfinder <- function(obj, sample_id, mpK){
  nCells <- ncol(obj)
  rate <- estimate_doublet_rate(nCells)
  nExp_poi <- round(rate * nCells)
  
  homotypic.prop <- modelHomotypic(obj$seurat_clusters)
  nExp_adj <- round(nExp_poi * (1 - homotypic.prop))
  
  obj <- doubletFinder(obj, PCs = PCs_use, pN = 0.25, pK = mpK, nExp = nExp_adj, sct = FALSE)
  
  df_col <- grep("^DF\\.classifications_", colnames(obj@meta.data), value = TRUE)
  df_col <- df_col[length(df_col)]
  obj$doublets <- obj@meta.data[[df_col]]
  
  # 可视化
  pdf(file.path(fig_dir, paste0(sample_id, "_umap_doublet.pdf")), width = 10, height = 4.5)
  print(
    DimPlot(obj, group.by = "doublets", shuffle = TRUE) + NoLegend() |
      DimPlot(obj, group.by = "doublets", split.by = "doublets", shuffle = TRUE) + NoLegend()
  )
  dev.off()
  
  list(
    obj = obj,
    stats = data.frame(
      sample_id = sample_id,
      nCells = nCells,
      doublet_rate = rate,
      nExp_poi = nExp_poi,
      homotypic_prop = round(homotypic.prop, 3),
      nExp_adj = nExp_adj,
      pK = mpK,
      n_doublet = sum(obj$doublets == "Doublet", na.rm = TRUE),
      n_singlet = sum(obj$doublets == "Singlet", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}

# -----------------------------
# main
# -----------------------------
rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop_msg("No .rds found in: ", in_dir)

sample_ids <- sub("\\.rds$", "", basename(rds_files))
rds_map <- setNames(rds_files, sample_ids)

# 选 pK 参考样本
if (pk_mode == "reuse_one") {
  if (is.null(pk_reference_sample)) {
    # 自动选“细胞数中位数”的样本
    sizes <- sapply(sample_ids, function(sid) {
      obj <- readRDS(rds_map[[sid]])
      ncol(obj)
    })
    pk_reference_sample <- names(sort(sizes))[ceiling(length(sizes)/2)]
  }
  message("[03] pK reference sample: ", pk_reference_sample)
}

# 先准备 pK（若 reuse_one）
mpK_global <- NA_real_
if (pk_mode == "reuse_one") {
  obj_ref <- readRDS(rds_map[[pk_reference_sample]])
  obj_ref <- prep_for_df(obj_ref)
  mpK_global <- find_best_pk(obj_ref, paste0(pk_reference_sample, "_REF"))
  message("[03] mpK_global = ", mpK_global)
  rm(obj_ref); gc()
}

all_stats <- list()

for (sid in sample_ids) {
  message("[03] ", sid)
  
  obj <- readRDS(rds_map[[sid]])
  obj <- prep_for_df(obj)
  
  mpK <- if (pk_mode == "per_sample") find_best_pk(obj, sid) else mpK_global
  
  res <- run_doubletfinder(obj, sid, mpK)
  obj2 <- res$obj
  all_stats[[sid]] <- res$stats
  
  # 保存：标注后的对象（或去除doublet）
  if (remove_doublets) {
    obj_out <- subset(obj2, subset = doublets == "Singlet")
  } else {
    obj_out <- obj2
  }
  saveRDS(obj_out, file.path(out_dir, paste0(sid, ".rds")))
  
  rm(obj, obj2, obj_out); gc()
}

stats_df <- bind_rows(all_stats)
write_csv(stats_df, file.path(tab_dir, "03_doublet_summary.csv"))
message("[OK] wrote: ", file.path(tab_dir, "03_doublet_summary.csv"))
message("[DONE] DoubletFinder finished.")
