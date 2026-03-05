#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE="$ROOT/data/raw/10x"
OUT="$ROOT/metadata/samples_auto.csv"

echo "sample_id,mrna_path,tcr_path,bcr_path,loom_path,has_mrna,has_tcr,has_bcr" > "$OUT"
for d in "$BASE"/*; do
  [[ -d "$d" ]] || continue
  sid="$(basename "$d")"
  mrna="$d/mrna"; tcr="$d/tcr"; bcr="$d/bcr"; loom="$d/loom"
  hm=0; ht=0; hb=0
  [[ -e "$mrna" ]] && hm=1
  [[ -e "$tcr"  ]] && ht=1
  [[ -e "$bcr"  ]] && hb=1
  echo "$sid,$mrna,$tcr,$bcr,$loom,$hm,$ht,$hb" >> "$OUT"
done
echo "[make_samples] wrote: $OUT"
