#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
})

message("[05] cluster cleanup: start")

# -----------------------------
# paths
# -----------------------------
in_rds  <- "data/processed/04_integrated/sce_harmony_integrated.rds"
out_dir <- "data/processed/05_cleaned"
fig_dir <- "results/figures/05_cleanup"
tab_dir <- "results/tables"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

stop_msg <- function(...) stop(paste0(...), call. = FALSE)

# -----------------------------
# USER: manual controls (edit here)
# -----------------------------

# 1) 手动剔除的 cluster（先跑一次看 QC dashboard，再填）
clusters_drop <- c(
  "8","14" #C8：doublets；C14：cycling T
  # example: "12","19"
)

# 2) 可选：细胞级阈值（不想用就都设为 NULL）
cell_filters <- list(
  nFeature_min = NULL,
  nFeature_max = NULL,
  nCount_min   = NULL,
  nCount_max   = NULL,
  mt_max       = NULL,
  decontX_max  = NULL
)

# 3) Doublet 差异化处理
remove_doublets_in_T_NK <- TRUE     # 只删 T/NK 的 DF doublet
remove_mixed_lineage_doublets <- TRUE  # 删“混合谱系 marker”的高置信 doublet（不分谱系）

# 混合判定阈值（在 log-normalized data 上）
mix_expr_cut <- 1.0   # 越大越严格；1.0 通常够用

# 4) SingleR 相关配置
run_singler <- TRUE
singler_mode <- "cluster"   # "cluster" 推荐；也可改 "cell" 但会很慢

# -----------------------------
# lineage scoring markers
# -----------------------------
sig_T      <- c("Trac","Cd3d","Cd3e","Lck","Cd247")
sig_NK     <- c("Ncr1","Klrd1","Klrk1","Nkg7")
sig_B      <- c("Cd79a","Cd79b","Cd19","Ms4a1","Cd74")
sig_Myelo  <- c("Lyz2","Tyrobp","Ctss","Itgam","Fcgr3","Csf1r")
sig_DC     <- c("Itgax","Zbtb46","Flt3","H2-Ab1","H2-Aa","Cd74")
sig_Plasma <- c("Jchain","Xbp1","Mzb1","Prdm1","Sdc1")
sig_Cycle  <- c("Mki67","Top2a","Stmn1","Pcna")

# 高置信“混合谱系”检测用 marker（少而硬）
mix_markers <- c("Trac","Cd79a","Ms4a1","Lyz2","Tyrobp","Itgax")

# -----------------------------
# helpers
# -----------------------------
safe_intersect <- function(x, genes) intersect(x, genes)

get_expr_df <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)
  
  # Seurat v5 优先 layer="data"，否则 fallback slot="data"
  mat <- NULL
  mat <- tryCatch({
    GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  }, error = function(e) NULL)
  
  if (is.null(mat)) {
    mat <- GetAssayData(obj, assay = "RNA", slot = "data")[genes, , drop = FALSE]
  }
  
  df <- as.data.frame(t(as.matrix(mat)))
  df$cell <- rownames(df)
  df
}

# -----------------------------
# 1) load integrated object
# -----------------------------
if (!file.exists(in_rds)) stop_msg("missing: ", in_rds)
sce <- readRDS(in_rds)

# ensure Idents
if (!"seurat_clusters" %in% colnames(sce@meta.data)) {
  # FindClusters 默认会有 seurat_clusters；如果没有，至少用 Idents
  sce$seurat_clusters <- as.character(Idents(sce))
}
Idents(sce) <- sce$seurat_clusters
sce$seurat_clusters <- as.character(sce$seurat_clusters)

# compute complexity metric (for QC)
sce$log10GenesPerUMI <- log10(sce$nFeature_RNA) / log10(sce$nCount_RNA)

# -----------------------------
# 2) add module scores for lineage diagnosis (fast)
# -----------------------------
genes_all <- rownames(sce)
sig_T_use      <- safe_intersect(sig_T, genes_all)
sig_NK_use     <- safe_intersect(sig_NK, genes_all)
sig_B_use      <- safe_intersect(sig_B, genes_all)
sig_Myelo_use  <- safe_intersect(sig_Myelo, genes_all)
sig_DC_use     <- safe_intersect(sig_DC, genes_all)
sig_Plasma_use <- safe_intersect(sig_Plasma, genes_all)
sig_Cycle_use  <- safe_intersect(sig_Cycle, genes_all)

# AddModuleScore 会新增列如 TScore1 等
sce <- AddModuleScore(sce, features = list(sig_T_use),      name = "TScore")
sce <- AddModuleScore(sce, features = list(sig_NK_use),     name = "NKScore")
sce <- AddModuleScore(sce, features = list(sig_B_use),      name = "BScore")
sce <- AddModuleScore(sce, features = list(sig_Myelo_use),  name = "MyeloScore")
sce <- AddModuleScore(sce, features = list(sig_DC_use),     name = "DCScore")
sce <- AddModuleScore(sce, features = list(sig_Plasma_use), name = "PlasmaScore")
sce <- AddModuleScore(sce, features = list(sig_Cycle_use),  name = "CycleScore")

# major_lineage：取 max score（粗分即可）
score_cols <- c("TScore1","NKScore1","BScore1","MyeloScore1","DCScore1","PlasmaScore1")
score_cols <- intersect(score_cols, colnames(sce@meta.data))

md <- sce@meta.data
md$major_lineage <- "Other"
if (length(score_cols) > 0) {
  score_mat <- as.matrix(md[, score_cols, drop = FALSE])
  max_idx <- max.col(score_mat, ties.method = "first")
  max_name <- colnames(score_mat)[max_idx]
  
  # 映射到更短标签
  map <- c(
    TScore1="T",
    NKScore1="NK",
    BScore1="B",
    MyeloScore1="Myeloid",
    DCScore1="DC",
    PlasmaScore1="Plasma"
  )
  md$major_lineage <- unname(map[max_name])
  
  # 低置信：如果 max score 很低，标记为 Other（避免乱判）
  max_val <- score_mat[cbind(seq_len(nrow(score_mat)), max_idx)]
  md$major_lineage[is.na(max_val) | max_val < 0.05] <- "Other"
}
sce@meta.data <- md
rm(md); gc()

# -----------------------------
# 3) SingleR for lineage diagnosis (cluster level)
# -----------------------------
if (run_singler) {
  message("[05] SingleR annotation (", singler_mode, ") ...")
  
  # Seurat -> SingleCellExperiment（SingleR 最稳的输入）
  sce_sce <- as.SingleCellExperiment(sce)
  
  ref_mr  <- celldex::MouseRNAseqData()
  ref_img <- celldex::ImmGenData()
  
  if (singler_mode == "cluster") {
    # 只对 cluster 做注释（推荐）
    pred <- SingleR(
      test = sce_sce,
      ref = list(MR = ref_mr, IMG = ref_img),
      labels = list(ref_mr$label.main, ref_img$label.main),
      clusters = sce$seurat_clusters
    )
    
    # pred 的行名就是 cluster id
    singler_cluster <- data.frame(
      seurat_clusters = rownames(pred),
      singleR_label = pred$labels,
      singleR_pruned = pred$pruned.labels,
      stringsAsFactors = FALSE
    )
    
    # 映射回每个细胞
    sce$singleR_cluster <- singler_cluster$singleR_label[match(sce$seurat_clusters, singler_cluster$seurat_clusters)]
    sce$singleR_cluster_pruned <- singler_cluster$singleR_pruned[match(sce$seurat_clusters, singler_cluster$seurat_clusters)]
    
    # 输出表
    readr::write_csv(singler_cluster, file.path(tab_dir, "05_singleR_cluster_labels.csv"))
    
  } else {
    # 逐细胞注释（非常慢）
    pred <- SingleR(
      test = sce_sce,
      ref = list(MR = ref_mr, IMG = ref_img),
      labels = list(ref_mr$label.main, ref_img$label.main)
    )
    sce$singleR_cell <- pred$labels
    sce$singleR_cell_pruned <- pred$pruned.labels
    readr::write_csv(
      data.frame(cell = colnames(sce), singleR = sce$singleR_cell, pruned = sce$singleR_cell_pruned),
      file.path(tab_dir, "05_singleR_cell_labels.csv")
    )
  }
  
  rm(sce_sce, ref_mr, ref_img, pred); gc()
}

# -----------------------------
# 4) cluster-level QC summary tables
# -----------------------------
meta <- sce@meta.data %>%
  rownames_to_column("cell")

# 基础 QC 指标列（存在就用）
qc_cols <- c("nFeature_RNA","nCount_RNA","mt_percent","decontX_contam","log10GenesPerUMI","ISGscore1","S.Score","G2M.Score")
qc_cols <- intersect(qc_cols, colnames(meta))

cluster_qc <- meta %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    # 中位数更稳
    across(all_of(qc_cols), ~ median(.x, na.rm = TRUE), .names = "median_{.col}"),
    doublet_frac = if ("doublets" %in% colnames(meta)) mean(doublets == "Doublet", na.rm = TRUE) else NA_real_,
    T_frac = mean(major_lineage == "T", na.rm = TRUE),
    NK_frac = mean(major_lineage == "NK", na.rm = TRUE),
    B_frac = mean(major_lineage == "B", na.rm = TRUE),
    Myeloid_frac = mean(major_lineage == "Myeloid", na.rm = TRUE),
    DC_frac = mean(major_lineage == "DC", na.rm = TRUE),
    Plasma_frac = mean(major_lineage == "Plasma", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write_csv(cluster_qc, file.path(tab_dir, "05_cluster_qc_summary.csv"))

# Phase 构成（如果有）
if ("Phase" %in% colnames(meta)) {
  cluster_phase <- meta %>%
    group_by(seurat_clusters, Phase) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()
  write_csv(cluster_phase, file.path(tab_dir, "05_cluster_phase_composition.csv"))
}

# sample / antigen / Group 分布（如果有）
if ("orig.ident" %in% colnames(meta)) {
  cluster_sample <- meta %>%
    group_by(seurat_clusters, orig.ident) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()
  write_csv(cluster_sample, file.path(tab_dir, "05_cluster_sample_composition.csv"))
}
if ("Group" %in% colnames(meta)) {
  cluster_group <- meta %>%
    group_by(seurat_clusters, Group) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()
  write_csv(cluster_group, file.path(tab_dir, "05_cluster_group_composition.csv"))
}

message("[05] wrote tables: 05_cluster_qc_summary.csv (+ phase/sample/group if available)")

# -----------------------------
# 5) QC dashboard plots (global)
# -----------------------------
# UMAP cluster
p_umap_cluster <- DimPlot(sce, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) +
  ggtitle("UMAP clusters (pre-clean)")
ggsave(file.path(fig_dir, "01_umap_clusters_preclean.png"), p_umap_cluster, width = 8, height = 6,dpi = 600)

# UMAP by lineage/doublet/phase
p_umap_lineage <- DimPlot(sce, reduction = "umap", group.by = "major_lineage", raster = TRUE) +
  ggtitle("UMAP by major_lineage")
ggsave(file.path(fig_dir, "02_umap_by_major_lineage.png"), p_umap_lineage, width = 8, height = 6)

if ("doublets" %in% colnames(sce@meta.data)) {
  p_umap_doublet <- DimPlot(sce, reduction = "umap", group.by = "doublets", raster = TRUE) +
    ggtitle("UMAP by DoubletFinder")
  ggsave(file.path(fig_dir, "03_umap_by_doublets.png"), p_umap_doublet, width = 8, height = 6)
}

if ("Phase" %in% colnames(sce@meta.data)) {
  p_umap_phase <- DimPlot(sce, reduction = "umap", group.by = "Phase", raster = TRUE) +
    ggtitle("UMAP by Phase")
  ggsave(file.path(fig_dir, "04_umap_by_phase.png"), p_umap_phase, width = 8, height = 6)
}

# QC feature plots
feat_show <- c("nFeature_RNA","nCount_RNA","mt_percent","decontX_contam","log10GenesPerUMI","ISGscore1")
feat_show <- intersect(feat_show, colnames(sce@meta.data))
if (length(feat_show) > 0) {
  p_feat <- FeaturePlot(sce, features = feat_show[1:min(6, length(feat_show))], ncol = 3, raster = TRUE)
  ggsave(file.path(fig_dir, "05_featureplot_qc_metrics.png"), p_feat, width = 10, height = 8)
}

# Violin plots by cluster (关键 QC 变量)
vln_feats <- c("nFeature_RNA","nCount_RNA","mt_percent","decontX_contam","log10GenesPerUMI","ISGscore1","CycleScore1")
vln_feats <- intersect(vln_feats, colnames(sce@meta.data))
if (length(vln_feats) > 0) {
  p_vln <- VlnPlot(sce, features = vln_feats, group.by = "seurat_clusters", pt.size = 0, ncol = 3)
  ggsave(file.path(fig_dir, "06_vlnplot_by_cluster.png"), p_vln, width = 20, height = 12)
}

# Cluster size + doublet fraction barplot
p_bar_n <- ggplot(cluster_qc, aes(x = reorder(seurat_clusters, n_cells), y = n_cells)) +
  geom_col() + coord_flip() + ggtitle("Cluster size (n_cells)")
ggsave(file.path(fig_dir, "07_cluster_size.png"), p_bar_n, width = 8, height = 10)

if (!all(is.na(cluster_qc$doublet_frac))) {
  p_bar_d <- ggplot(cluster_qc, aes(x = reorder(seurat_clusters, doublet_frac), y = doublet_frac)) +
    geom_col() + coord_flip() + ggtitle("Doublet fraction per cluster")
  ggsave(file.path(fig_dir, "08_cluster_doublet_fraction.png"), p_bar_d, width = 8, height = 10)
}

# DotPlot: lineage + contamination markers by cluster (帮助人工判 cluster)
dot_markers <- list(
  T = c("Trac","Cd3d","Cd3e"),
  B = c("Cd79a","Ms4a1","Cd74"),
  Myeloid = c("Lyz2","Tyrobp","Ctss"),
  DC = c("Itgax","Flt3","Zbtb46"),
  Plasma = c("Jchain","Xbp1","Mzb1"),
  Cycle = c("Mki67","Top2a"),
  RBC_Plt = c("Hbb-bs","Hba-a1","Pf4","Ppbp")
)
dot_genes <- unique(unlist(dot_markers))
dot_genes <- intersect(dot_genes, rownames(sce))
if (length(dot_genes) > 0) {
  p_dot <- DotPlot(sce, features = dot_markers, group.by = "seurat_clusters") +
    RotatedAxis() + theme_bw() + theme(panel.grid = element_blank())
  ggsave(file.path(fig_dir, "09_dotplot_lineage_markers_by_cluster.png"), p_dot, width = 14, height = 8)
}

# -----------------------------
# 6) build "mixed lineage" flag (cell-level)
# -----------------------------
mix_df <- get_expr_df(sce, mix_markers)
if (!is.null(mix_df)) {
  # 判断：Trac 与 (Lyz2/Tyrobp) 或 (Cd79a/Ms4a1) 同时高；或 B 与 Myeloid 同时高
  has <- function(g) if (g %in% colnames(mix_df)) mix_df[[g]] else 0
  
  trac  <- has("Trac")
  cd79a <- has("Cd79a")
  ms4a1 <- has("Ms4a1")
  lyz2  <- has("Lyz2")
  tyrobp<- has("Tyrobp")
  itgax <- has("Itgax")
  
  mix_flag <- (trac > mix_expr_cut & (lyz2 > mix_expr_cut | tyrobp > mix_expr_cut | itgax > mix_expr_cut)) |
    (trac > mix_expr_cut & (cd79a > mix_expr_cut | ms4a1 > mix_expr_cut)) |
    ((cd79a > mix_expr_cut | ms4a1 > mix_expr_cut) & (lyz2 > mix_expr_cut | tyrobp > mix_expr_cut))
  names(mix_flag) <- mix_df$cell
  sce$mix_lineage_flag <- FALSE
  sce$mix_lineage_flag[names(mix_flag)] <- mix_flag
} else {
  sce$mix_lineage_flag <- FALSE
}

# -----------------------------
# 7) apply cleanup rules
# -----------------------------
meta2 <- sce@meta.data %>% rownames_to_column("cell")

keep <- rep(TRUE, nrow(meta2))
reason <- rep("", nrow(meta2))

# (A) cluster drop
if (length(clusters_drop) > 0) {
  drop_idx <- meta2$seurat_clusters %in% clusters_drop
  keep[drop_idx] <- FALSE
  reason[drop_idx] <- paste0(reason[drop_idx], "drop_cluster;")
}

# (B) cell-level numeric filters
apply_num_filter <- function(col, op, thr, tag) {
  if (is.null(thr)) return(invisible(NULL))
  if (!col %in% colnames(meta2)) return(invisible(NULL))
  x <- meta2[[col]]
  idx <- op(x, thr)
  keep[idx] <<- FALSE
  reason[idx] <<- paste0(reason[idx], tag, ";")
}

apply_num_filter("nFeature_RNA", `<`, cell_filters$nFeature_min, "nFeature_low")
apply_num_filter("nFeature_RNA", `>`, cell_filters$nFeature_max, "nFeature_high")
apply_num_filter("nCount_RNA",   `<`, cell_filters$nCount_min,   "nCount_low")
apply_num_filter("nCount_RNA",   `>`, cell_filters$nCount_max,   "nCount_high")
apply_num_filter("mt_percent",   `>`, cell_filters$mt_max,       "mt_high")
apply_num_filter("decontX_contam", `>`, cell_filters$decontX_max, "decontX_high")

# (C) remove DF doublets in T/NK only
if (remove_doublets_in_T_NK && ("doublets" %in% colnames(meta2))) {
  idx <- (meta2$doublets == "Doublet") & (meta2$major_lineage %in% c("T","NK"))
  keep[idx] <- FALSE
  reason[idx] <- paste0(reason[idx], "DF_doublet_TNK;")
}

# (D) remove mixed-lineage (high confidence) doublets (or even all mixed cells)
if (remove_mixed_lineage_doublets && ("mix_lineage_flag" %in% colnames(meta2))) {
  idx <- meta2$mix_lineage_flag
  keep[idx] <- FALSE
  reason[idx] <- paste0(reason[idx], "mixed_lineage;")
}

cells_keep <- meta2$cell[keep]
cells_drop <- meta2$cell[!keep]

# 输出被删除细胞表（方便你回溯）
removed_df <- meta2[!keep, c("cell","orig.ident","batch","age_week","antigen","sex","Group",
                             "seurat_clusters","major_lineage","doublets",
                             "nFeature_RNA","nCount_RNA","mt_percent","decontX_contam","log10GenesPerUMI",
                             "ISGscore1","S.Score","G2M.Score","Phase","mix_lineage_flag"), drop = FALSE]
removed_df$reason <- reason[!keep]

write_csv(removed_df, file.path(tab_dir, "05_cells_removed.csv"))

# 统计删除原因
reason_summary <- data.frame(reason = removed_df$reason) %>%
  mutate(reason = ifelse(is.na(reason) | reason == "", "unknown", reason)) %>%
  count(reason, sort = TRUE)
write_csv(reason_summary, file.path(tab_dir, "05_cells_removed_reason_summary.csv"))

message("[05] keep cells: ", length(cells_keep), " / drop cells: ", length(cells_drop))

# 子集得到 cleaned 对象
sce_clean <- subset(sce, cells = cells_keep)

# 保存清洗后对象
saveRDS(sce_clean, file.path(out_dir, "sce_cleaned.rds"))

# 同时输出清洗后基本对照图
p_post1 <- DimPlot(sce_clean, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) +
  ggtitle("UMAP clusters (post-clean)")
ggsave(file.path(fig_dir, "10_umap_clusters_postclean.png"), p_post1, width = 9, height = 6)

p_post2 <- DimPlot(sce_clean, reduction = "umap", group.by = "major_lineage", raster = TRUE) +
  ggtitle("UMAP by major_lineage (post-clean)")
ggsave(file.path(fig_dir, "11_umap_lineage_postclean.png"), p_post2, width = 8, height = 6)

if ("doublets" %in% colnames(sce_clean@meta.data)) {
  p_post3 <- DimPlot(sce_clean, reduction = "umap", group.by = "doublets", raster = TRUE) +
    ggtitle("UMAP by DoubletFinder (post-clean)")
  ggsave(file.path(fig_dir, "12_umap_doublets_postclean.png"), p_post3, width = 7, height = 5)
}

message("[OK] saved cleaned object: ", file.path(out_dir, "sce_cleaned.rds"))
message("[05] cluster cleanup: done")