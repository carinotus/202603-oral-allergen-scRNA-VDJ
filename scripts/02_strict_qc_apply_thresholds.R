#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# user paths
# -----------------------------
in_dir  <- "data/processed/01_preQC"
out_dir <- "data/processed/02_QC"
fig_dir <- "results/figures/02_strict_qc"
tab_dir <- "results/tables"
thr_fp  <- "metadata/qc_thresholds.csv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

stop_msg <- function(...) stop(paste0(...), call. = FALSE)

# -----------------------------
# helpers
# -----------------------------
q2 <- function(x, probs = c(0.03, 0.97)) {
  if (is.null(x) || length(x) == 0) return(c(NA_real_, NA_real_))
  if (all(is.na(x))) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

calc_qc_quantiles <- function(obj) {
  md <- obj@meta.data
  getv <- function(k) if (k %in% colnames(md)) md[[k]] else rep(NA_real_, nrow(md))
  
  q_nf <- q2(getv("nFeature_RNA"))
  q_nc <- q2(getv("nCount_RNA"))
  q_mt <- q2(getv("mt_percent"))
  q_dx <- q2(getv("decontX_contam"))
  
  data.frame(
    n_before = ncol(obj),
    nFeature_q03 = q_nf[1], nFeature_q97 = q_nf[2],
    nCount_q03   = q_nc[1], nCount_q97   = q_nc[2],
    mt_q03       = q_mt[1], mt_q97       = q_mt[2],
    decontX_q03  = q_dx[1], decontX_q97  = q_dx[2],
    stringsAsFactors = FALSE
  )
}

plot_threshold_scatter <- function(obj, cut, title_prefix = "") {
  f_low   <- as.numeric(cut$f_low)
  f_high  <- as.numeric(cut$f_high)
  c_low   <- as.numeric(cut$c_low)
  c_high  <- as.numeric(cut$c_high)
  mt_high <- as.numeric(cut$mt_high)
  
  p11 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept = f_low,  linetype = "dashed") +
    geom_hline(yintercept = f_high) +
    geom_vline(xintercept = c_low,  linetype = "dashed") +
    geom_vline(xintercept = c_high) +
    ggtitle(paste0(title_prefix, "  nCount vs nFeature"))
  
  p12 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "mt_percent") +
    geom_hline(yintercept = mt_high) +
    geom_vline(xintercept = c_low,  linetype = "dashed") +
    geom_vline(xintercept = c_high) +
    ggtitle(paste0(title_prefix, "  nCount vs mt%"))
  
  p11 + p12 + plot_layout(ncol = 1)
}

apply_strict_qc <- function(obj, cut) {
  # 强制数值化（避免阈值读成字符）
  f_low   <- as.numeric(cut$f_low)
  f_high  <- as.numeric(cut$f_high)
  c_low   <- as.numeric(cut$c_low)
  c_high  <- as.numeric(cut$c_high)
  mt_high <- as.numeric(cut$mt_high)
  
  subset(
    obj,
    subset =
      nFeature_RNA > f_low  &
      nFeature_RNA < f_high &
      nCount_RNA   > c_low  &
      nCount_RNA   < c_high &
      mt_percent   < mt_high
  )
}

# -----------------------------
# 0) index rds + load cache once
# -----------------------------
rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop_msg("No .rds found in: ", in_dir)

sample_ids <- sub("\\.rds$", "", basename(rds_files))
rds_map <- setNames(rds_files, sample_ids)

message("[02] loading all preQC objects into memory (seurat_cache)...")
seurat_cache <- setNames(vector("list", length(sample_ids)), sample_ids)
for (sid in sample_ids) {
  message("[02] load: ", sid)
  seurat_cache[[sid]] <- readRDS(rds_map[[sid]])
}
message("[02] cache ready: ", length(seurat_cache), " objects")

# -----------------------------
# 1) quantile table (for decision)
# -----------------------------
message("[02] generating quantile table (preQC)...")
quant_rows <- lapply(names(seurat_cache), function(sid) {
  obj <- seurat_cache[[sid]]
  qdf <- calc_qc_quantiles(obj)
  qdf$sample_id <- sid
  qdf
})
quant_df <- dplyr::bind_rows(quant_rows)
quant_df <- quant_df[, c("sample_id", setdiff(colnames(quant_df), "sample_id"))]

write_csv(quant_df, file.path(tab_dir, "02_qc_quantiles_preQC.csv"))
message("[OK] wrote: ", file.path(tab_dir, "02_qc_quantiles_preQC.csv"))

# -----------------------------
# 2) read thresholds
# -----------------------------
message("[02] reading thresholds: ", thr_fp)
thr <- read_csv(thr_fp, show_col_types = FALSE, comment = "#")
need_cols <- c("sample_id","f_low","f_high","c_low","c_high","mt_high")
if (!all(need_cols %in% colnames(thr))) {
  stop_msg("qc_thresholds.csv must contain columns: ",
           paste(need_cols, collapse = ", "))
}
thr <- thr %>% mutate(sample_id = as.character(sample_id))

thr_bad <- thr %>% filter(if_any(all_of(need_cols[-1]), ~ is.na(.x)))
if (nrow(thr_bad) > 0) {
  message("ERROR: some thresholds are NA, please fill them. Examples:")
  print(head(thr_bad, 10))
  stop_msg("Threshold table contains NA values.")
}

missing_obj <- setdiff(thr$sample_id, names(seurat_cache))
if (length(missing_obj) > 0) {
  stop_msg("These sample_id are in qc_thresholds.csv but missing in seurat_cache:\n  ",
           paste(missing_obj, collapse = ", "))
}

# -----------------------------
# 3) strict QC using cache + save
# -----------------------------
# QC scatter plot
for (i in seq_len(nrow(thr))) {
  sid <- thr$sample_id[i]
  cut <- thr[i, ]
  
  message("[02] QC scatter plot: ", sid)
  
  obj <- seurat_cache[[sid]]
  n_before <- ncol(obj)
  
  # plot thresholds
  p <- plot_threshold_scatter(obj, cut, title_prefix = sid)
  pdf(file.path(fig_dir, paste0(sid, "_thresholds.pdf")), width = 7.5, height = 9)
  print(p)
  dev.off()
}

# strict QC
qc_summary <- vector("list", nrow(thr))
names(qc_summary) <- thr$sample_id

for (i in seq_len(nrow(thr))) {
  sid <- thr$sample_id[i]
  cut <- thr[i, ]
  
  message("[02] strict QC: ", sid)
  
  obj <- seurat_cache[[sid]]
  n_before <- ncol(obj)
  
  # strict QC
  obj2 <- apply_strict_qc(obj, cut)
  n_after <- ncol(obj2)
  
  saveRDS(obj2, file.path(out_dir, paste0(sid, ".rds")))
  
  qc_summary[[sid]] <- data.frame(
    sample_id = sid,
    n_before = n_before,
    n_after = n_after,
    kept_frac = ifelse(n_before > 0, n_after / n_before, NA_real_),
    f_low = cut$f_low, f_high = cut$f_high,
    c_low = cut$c_low, c_high = cut$c_high,
    mt_high = cut$mt_high,
    stringsAsFactors = FALSE
  )
}

qc_df <- dplyr::bind_rows(qc_summary)
write_csv(qc_df, file.path(tab_dir, "02_qc_summary.csv"))
message("[OK] wrote: ", file.path(tab_dir, "02_qc_summary.csv"))
message("[DONE] strict QC finished.")
