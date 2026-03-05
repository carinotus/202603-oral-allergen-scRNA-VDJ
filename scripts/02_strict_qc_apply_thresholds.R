#!/usr/bin/env Rscript
message("[02] soup decontam: start")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

thr <- read.csv("metadata/qc_thresholds.csv", stringsAsFactors = FALSE, comment.char = "#")
stopifnot(all(c("sample_id","f_low","f_high","c_low","c_high","mt_high") %in% colnames(thr)))

dir.create("data/processed/02_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/02_qc_after", recursive = TRUE, showWarnings = FALSE)

apply_one <- function(sample_id, cut){
  obj <- readRDS(file.path("data/processed/01_preQC", paste0(sample_id, ".rds")))
  
  obj2 <- subset(obj, subset =
                   nFeature_RNA > cut$f_low &
                   nFeature_RNA < cut$f_high &
                   nCount_RNA   > cut$c_low &
                   nCount_RNA   < cut$c_high &
                   mt_percent   < cut$mt_high)
  
  saveRDS(obj2, file.path("data/processed/02_QC", paste0(sample_id, ".rds")))
  invisible(obj2)
}

for(i in seq_len(nrow(thr))){
  sid <- thr$sample_id[i]
  cut <- thr[i,]
  message("[02] strict QC: ", sid)
  apply_one(sid, cut)
}