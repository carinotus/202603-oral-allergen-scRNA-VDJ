#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE="$ROOT/data/raw/10x"
OUT="$ROOT/metadata/samples.csv"

# header
echo "sample_id,batch,age_week,antigen,sex,replicate,mrna_type,mrna_path,tcr_path,bcr_path,loom_path,notes" > "$OUT"

# helper: infer batch/age/antigen/sex from sample_id pattern
infer_meta(){
  local sid="$1"
  local batch=""; local age=""; local antigen=""; local sex=""; local rep=""; local notes=""

  # ---- heuristic rules (你后面可手动改) ----
  # 3W_PBS_HDM: CONT_*, OR_*  (按你目前目录推断)
  if [[ "$sid" =~ ^CONT_ || "$sid" =~ ^OR_ ]]; then
    batch="3W_PBS_HDM"
    age="3"
    antigen="HDM"   # 如果你这里其实是 PBS vs HDM 的混合，请后面手动改
  elif [[ "$sid" =~ ^(PBS|OVA)_(F|M)_ ]]; then
    batch="3W_PBS_OVA"
    age="3"
    antigen="${sid%%_*}"          # PBS 或 OVA
    sex="$(echo "$sid" | cut -d'_' -f2)"  # F / M
  elif [[ "$sid" =~ ^LC_ ]]; then
    batch="5W_PBS_HDM_OVA"
    age="5"
    antigen="UNK"   # 这个批次你需要自己补：PBS/HDM/OVA
  else
    batch="UNK"; age=""; antigen=""; notes="check_id_pattern"
  fi

  echo "$batch,$age,$antigen,$sex,$rep,$notes"
}

# loop samples
for d in "$BASE"/*; do
  [[ -d "$d" ]] || continue
  sid="$(basename "$d")"

  mrna="$d/mrna"
  tcr="$d/tcr"
  bcr="$d/bcr"
  loom="$d/loom"

  # detect mrna type
  mrna_type=""
  if [[ -d "$mrna" ]]; then
    # mtx folder: matrix.mtx(.gz) + barcodes + features
    if ls "$mrna"/matrix.mtx* 1>/dev/null 2>&1; then
      mrna_type="mtx_dir"
    else
      mrna_type="dir"
    fi
  elif [[ -f "$d/mrna.h5" ]]; then
    mrna="$d/mrna.h5"
    mrna_type="h5"
  else
    mrna_type="missing"
  fi

  meta="$(infer_meta "$sid")"
  IFS=',' read -r batch age antigen sex rep notes <<< "$meta"

  echo "$sid,$batch,$age,$antigen,$sex,$rep,$mrna_type,$mrna,$tcr,$bcr,$loom,$notes" >> "$OUT"
done

echo "[make_metadata] wrote: $OUT"
