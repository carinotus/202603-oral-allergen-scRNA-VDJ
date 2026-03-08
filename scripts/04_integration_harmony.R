#!/usr/bin/env Rscript
message("[04] integration: start")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(harmony)
})

# -----------------------------
# paths
# -----------------------------
in_dir  <- "data/processed/03_doublet"        # 不删doublet就用这里
out_dir <- "data/processed/04_integrated"
fig_dir <- "results/figures/04_integration_harmony"
tab_dir <- "results/tables"
samples_fp <- "metadata/samples.csv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# user settings 
# -----------------------------
npcs <- 50
dims_use <- 1:30
cluster_res <- 0.8

# Harmony batch variable
batch_var <- "orig.ident"   # sample-level batch correction

# regress confounders BEFORE PCA (optional)
# 先只回归 decontX_contam,ISG/CellCycle 先不回归？先画图判断
regress_decontX <- FALSE
regress_mt      <- FALSE
regress_cc      <- FALSE     # S.Score/G2M.Score
regress_isg     <- TRUE     # ISG module score

# ISG gene set (mouse)
isg_genes <- c("Irf7","Stat1","Stat2","Irf9","Isg15","Ifit1","Ifit2","Ifit3","Mx1","Oasl1","Oas1a","Oas2",
               "Usp18","Rsad2","Ifi35","Ifi44","Ifi47","Ddx58","Ifih1","Zbp1","Gbp2","Gbp3","Gbp4","Gbp5",
               "Bst2","Rtp4","Parp9","Parp14","Pml")

# quick marker sanity check
apc_markers <- c("Lyz2","Tyrobp","Ctss","Itgam","Itgax","Cd74","H2-Ab1","H2-Aa","Flt3","Zbtb46","Adgre1","Csf1r")

# -----------------------------
# helpers
# -----------------------------
stop_msg <- function(...) stop(paste0(...), call. = FALSE)

# -----------------------------
# 1) load all objects
# -----------------------------
rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop_msg("No .rds found in: ", in_dir)

sample_ids <- sub("\\.rds$", "", basename(rds_files))
obj_list <- setNames(vector("list", length(sample_ids)), sample_ids)

message("[04] loading objects ...")
for (sid in sample_ids) {
  message("[04] load: ", sid)
  obj <- readRDS(file.path(in_dir, paste0(sid, ".rds")))
  
  # 强制统一 sample metadata
  obj$orig.ident <- sid
  obj$samples <- sid
  
  obj_list[[sid]] <- obj
}

# -----------------------------
# 2) merge + JoinLayers + add metadata
# -----------------------------
message("[04] merging ...")
sce <- Reduce(function(x, y) merge(x = x, y = y, merge.data = FALSE), obj_list)
colnames(sce@meta.data)
saveRDS(sce, file.path(out_dir, paste0("sce_before_joinlaters", ".rds")))
sce <- JoinLayers(sce)

# 清理 meta（只保留关心的列）
keep_cols <- c(
  "orig.ident","samples",
  "nCount_RNA","nFeature_RNA","mt_percent","ribo_percent","decontX_contam",
  "doublets"
)
keep_cols <- intersect(keep_cols, colnames(sce@meta.data))
sce@meta.data <- sce@meta.data[, keep_cols, drop = FALSE]


# add metadata from files
samples_fp <- "metadata/samples.csv"
samples_meta <- readr::read_csv(samples_fp, show_col_types = FALSE)

# 只取需要写进每个细胞的列（路径列不要 join 到每个细胞，太大）
samples_meta_small <- samples_meta %>%
  dplyr::select(sample_id, batch, age_week, antigen, sex, replicate, notes)

# sample_id 必须唯一（否则 join 会把细胞行数翻倍）
if (anyDuplicated(samples_meta_small$sample_id)) {
  dup <- unique(samples_meta_small$sample_id[duplicated(samples_meta_small$sample_id)])
  stop("metadata/samples.csv sample_id not unique: ", paste(dup, collapse = ", "))
}

# 把每个细胞的 meta 拿出来，按 sample_id（=orig.ident）做 join，再放回去
md <- sce@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(sample_id = .data$orig.ident) %>%
  dplyr::left_join(samples_meta_small, by = "sample_id") %>%
  tibble::column_to_rownames("cell")

# 检查有没有没匹配到的样本（NA）
if (any(is.na(md$batch))) {
  unmapped <- unique(md$sample_id[is.na(md$batch)])
  warning("These sample_id are missing in metadata/samples.csv: ",
          paste(unmapped, collapse = ", "))
}

# 类型整理（可选但推荐）
md$batch   <- factor(md$batch)
md$age_week <- as.integer(md$age_week)
md$antigen <- factor(md$antigen, levels = c("PBS","HDM","OVA"))
md$sex     <- factor(md$sex, levels = c("F","M","Unknown"))

# 生成一个常用分组列（可选）：age+antigen
md$Group <- paste0(md$age_week, "W_", as.character(md$antigen))
md$Group <- factor(md$Group)
colnames(md)
table(md$orig.ident,md$Group)

sce@meta.data <- md
rm(md, samples_meta, samples_meta_small); gc()
view(sce@meta.data)


# 保存样本细胞数统计
cell_count <- as.data.frame(table(sce$orig.ident))
colnames(cell_count) <- c("sample_id","n_cells")
write_csv(cell_count, file.path(tab_dir, "04_cells_per_sample_pre_integration.csv"))

p_bar <- ggplot(cell_count, aes(x = reorder(sample_id, n_cells), y = n_cells)) +
  geom_col() + coord_flip() + ggtitle("Cells per sample (pre-integration)")
ggsave(file.path(fig_dir, "00_cells_per_sample.png"), p_bar, width = 8, height = 6)

# -----------------------------
# 3) optional: cell cycle + ISG score (for diagnosis / optional regression)
# -----------------------------
# Cell cycle scoring（可用于诊断；是否回归由开关控制）
message("[04] computing CellCycle scores (for diagnosis) ...")
s.genes = alias2SymbolTable(str_to_title(cc.genes$s.genes), species="Mm")
g2m.genes = alias2SymbolTable(str_to_title(cc.genes$g2m.genes), species="Mm")

sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# ISG module score（可用于诊断；是否回归由开关控制）
isg_genes_use <- intersect(isg_genes, rownames(sce))
if (length(isg_genes_use) >= 5) {
  sce <- AddModuleScore(sce, features = list(isg_genes_use), name = "ISGscore")
} else {
  message("[04] WARN: too few ISG genes found in object, skip ISGscore.")
}

# -----------------------------
# 4) Normalize/HVG/Scale/PCA (regress optional)
# -----------------------------
message("[04] Normalize/HVG/Scale/PCA ...")
sce <- NormalizeData(sce, verbose = FALSE)
sce <- FindVariableFeatures(sce, nfeatures = 2000, verbose = FALSE)

vars_regress <- c()
if (regress_decontX && "decontX_contam" %in% colnames(sce@meta.data)) vars_regress <- c(vars_regress, "decontX_contam")
if (regress_mt      && "mt_percent"     %in% colnames(sce@meta.data)) vars_regress <- c(vars_regress, "mt_percent")
if (regress_cc      && all(c("S.Score","G2M.Score") %in% colnames(sce@meta.data))) vars_regress <- c(vars_regress, "S.Score","G2M.Score")
if (regress_isg     && "ISGscore1"      %in% colnames(sce@meta.data)) vars_regress <- c(vars_regress, "ISGscore1")

message("[04] vars.to.regress = ", ifelse(length(vars_regress)==0, "NONE", paste(vars_regress, collapse = ", ")))
sce <- ScaleData(sce, vars.to.regress = vars_regress, verbose = FALSE)
sce <- RunPCA(sce, npcs = npcs, verbose = FALSE)

ggsave(file.path(fig_dir, "01_elbowplot.png"), ElbowPlot(sce, ndims = npcs), width = 6, height = 4)
ggsave(file.path(fig_dir, "02_pca_by_sample.png"),
       DimPlot(sce, reduction = "pca", group.by = "orig.ident", raster = TRUE) + NoLegend(),
       width = 8, height = 6)

# -----------------------------
# 5) Harmony (batch correction on PCA)
# -----------------------------
if (!batch_var %in% colnames(sce@meta.data)) stop_msg("batch var not found: ", batch_var)

message("[04] RunHarmony by ", batch_var, " ...")
sce <- RunHarmony(
  object = sce,
  group.by.vars = batch_var,
  # reduction = "pca",
  # dims.use = dims_use,
  # assay.use = "RNA",
  # verbose = TRUE
)
ggsave(file.path(fig_dir, "02_harmony_by_sample.png"),
       DimPlot(sce, reduction = "harmony", group.by = "orig.ident", raster = TRUE) + NoLegend(),
       width = 8, height = 6)

# -----------------------------
# 6) neighbors/clusters/umap on harmony
# -----------------------------
message("[04] clustering/umap on harmony ...")

sce <- FindNeighbors(sce, reduction = "harmony", dims = dims_use, verbose = FALSE)
sce <- FindClusters(sce, resolution = cluster_res, verbose = FALSE)
sce <- RunUMAP(sce, reduction = "harmony", dims = dims_use, verbose = FALSE) #此刻拜一拜

# -----------------------------
# 7) sanity-check plots
# -----------------------------
ggsave(file.path(fig_dir, "03_umap_by_sample.png"),
       DimPlot(sce, reduction = "umap", group.by = "orig.ident", raster = TRUE) + NoLegend() +
         ggtitle("UMAP (Harmony) by sample"),
       width = 9, height = 6)

pdf(file.path(fig_dir,"03_umap_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(sce, reduction = "umap", group.by = "orig.ident", raster = TRUE) +
        ggtitle("UMAP (Harmony) by sample")
      )
dev.off()

if ("doublets" %in% colnames(sce@meta.data)) {
  ggsave(file.path(fig_dir, "04_umap_by_doublets.png"),
         DimPlot(sce, reduction = "umap", group.by = "doublets", raster = TRUE) +
           ggtitle("UMAP (Harmony) by DoubletFinder"),
         width = 7, height = 5)
  pdf(file.path(fig_dir,"04_umap_by_doublets.pdf"), width = 8, height = 6)
  print(DimPlot(sce, reduction = "umap", group.by = "doublets", raster = TRUE) +
          ggtitle("UMAP (Harmony) by DoubletFinder")
  )
  dev.off()
}

ggsave(file.path(fig_dir, "05_umap_clusters.png"),
       DimPlot(sce, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) +
         ggtitle(paste0("UMAP clusters res=", cluster_res)),
       width = 9, height = 6)
pdf(file.path(fig_dir,"05_umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(sce, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) +
        ggtitle(paste0("UMAP clusters res=", cluster_res))
)
dev.off()

# Cell cycle / ISG diagnosis plots (do not imply regression unless toggled)
ggsave(file.path(fig_dir, "06_umap_by_phase.png"),
       DimPlot(sce, reduction = "umap", group.by = "Phase", raster = TRUE) +
         ggtitle("UMAP by cell cycle Phase"),
       width = 7, height = 5)
pdf(file.path(fig_dir,"06_umap_by_phase.pdf"),width = 8,height = 6)
print(DimPlot(sce,reduction = "umap",group.by = "Phase",raster = TRUE) + 
        ggtitle("UMAP by cell cycle Phase"))
dev.off()

if ("ISGscore1" %in% colnames(sce@meta.data)) {
  ggsave(file.path(fig_dir, "07_feature_isgscore.png"),
         FeaturePlot(sce, features = "ISGscore1", raster = TRUE) +
           ggtitle("ISG module score (diagnostic)"),
         width = 7, height = 5)
}
pdf(file.path(fig_dir,"07_feature_isgscore.pdf"),width = 8,height = 6)
print(FeaturePlot(sce, features = "ISGscore1", raster = TRUE) +
        ggtitle("ISG module score (diagnostic)")
      )
dev.off()

# APC markers quick check
apc_use <- intersect(apc_markers, rownames(sce))
if (length(apc_use) > 0) {
  ggsave(file.path(fig_dir, "08_featureplot_apc_markers.png"),
         FeaturePlot(sce, features = apc_use[1:min(9, length(apc_use))], ncol = 3, raster = TRUE),
         width = 10, height = 9)
}

# check HDM VS OVA
pdf(file.path(fig_dir,"09_UMAP spilt by antigen.pdf"),width = 20, height = 6)
print(DimPlot(sce, split.by = "antigen", raster = TRUE, label = T) + ggtitle("UMAP clusters spilt by antigen"))
dev.off()



# -----------------------------
# 8) save
# -----------------------------
saveRDS(sce, file.path(out_dir, "sce_harmony_integrated.rds"))
message("[OK] saved: ", file.path(out_dir, "sce_harmony_integrated.rds"))
message("[DONE] integration finished.")