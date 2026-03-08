#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN_DIR="$ROOT/data/processed/01_preQC"
OUT="$ROOT/metadata/qc_thresholds.csv"

mkdir -p "$ROOT/metadata"

# header
echo "sample_id,f_low,f_high,c_low,c_high,mt_high" > "$OUT"

# 从 01_preQC 的 rds 文件名自动提取 sample_id
# 你可以先填默认值（留空也行），这里我给留空，完全由你手填
for f in "$IN_DIR"/*.rds; do
  sid="$(basename "$f" .rds)"
  echo "$sid,,,,," >> "$OUT"
done

echo "[ok] wrote template: $OUT"
echo "edit it, then run: Rscript scripts/02_strict_qc_apply_thresholds.R"
