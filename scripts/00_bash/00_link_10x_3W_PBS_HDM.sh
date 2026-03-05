#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST_BASE="$ROOT/data/raw/10x"

MRNA_BASE="$ROOT/data/raw/3W_PBS_HDM/mrnaseq/unpacked/Summary/1_Cellranger_result"
TCR_BASE="$ROOT/data/raw/3W_PBS_HDM/tcrseq/unpacked/Summary/1_Cellranger_result"
BCR_BASE="$ROOT/data/raw/3W_PBS_HDM/bcrseq/unpacked/Summary/1_Cellranger_result"
LOOM_BASE="$ROOT/data/raw/3W_PBS_HDM/loom_file"

say(){ echo -e "[3W_HDM_link] $*"; }
safe_ln(){
  local src="$1" dst="$2"
  [[ -e "$dst" ]] && return 0
  ln -s "$src" "$dst"
}

[[ -d "$MRNA_BASE" ]] || { say "ERROR: missing $MRNA_BASE (did you extract mrnaseq Summary.tar.gz ?)"; exit 1; }
mkdir -p "$DEST_BASE"

# 以 mRNA 的样本目录为主，逐个建统一入口
while IFS= read -r -d '' sample_dir; do
  sid="$(basename "$sample_dir")"
  out="$DEST_BASE/$sid"
  mkdir -p "$out"

  # mRNA 矩阵：filtered_feature_bc_matrix 目录（里面是 barcodes/features/matrix）
  if [[ -d "$sample_dir/filtered_feature_bc_matrix" ]]; then
    safe_ln "$sample_dir/filtered_feature_bc_matrix" "$out/mrna"
  else
    say "WARN: $sid missing filtered_feature_bc_matrix/"
  fi

  # TCR/BCR：把整个 sample 文件夹链接过来（里面有 contig/统计/web_summary 等）
  [[ -d "$TCR_BASE/$sid" ]] && safe_ln "$TCR_BASE/$sid" "$out/tcr"
  [[ -d "$BCR_BASE/$sid" ]] && safe_ln "$BCR_BASE/$sid" "$out/bcr"

  # loom（如果你后面想用 loom 快速预览/导入）
  [[ -d "$LOOM_BASE/$sid" ]] && safe_ln "$LOOM_BASE/$sid" "$out/loom"

done < <(find "$MRNA_BASE" -mindepth 1 -maxdepth 1 -type d -print0)

say "Done. Example:"
ls -lah "$DEST_BASE" | head
