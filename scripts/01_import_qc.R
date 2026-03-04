#!/usr/bin/env Rscript
# Step 01: import + QC (per batch or per sample)
# Input: metadata/samples.csv + 10x matrices
# Output: QC plots + milestone Seurat objects in data/processed/

message("[01] import + QC: start")
