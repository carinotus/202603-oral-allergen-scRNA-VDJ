#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_DIRS=(
  "$ROOT/data/raw/3W_PBS_OVA/unpacked/Summary"
  "$ROOT/data/raw/5W_PBS_HDM_OVA/unpacked/Summary"
)
DEST_BASE="$ROOT/data/raw/10x"

say(){ echo -e "[link10x] $*"; }
safe_ln(){
  local src="$1" dst="$2"
  [[ -e "$dst" ]] && return 0
  ln -s "$src" "$dst"
}

mkdir -p "$DEST_BASE"

for summary in "${SRC_DIRS[@]}"; do
  [[ -d "$summary" ]] || { say "WARN: missing $summary (skip)"; continue; }

  while IFS= read -r -d '' sample_dir; do
    sample_id="$(basename "$sample_dir")"
    dest="$DEST_BASE/$sample_id"
    mkdir -p "$dest"

    # mRNA
    [[ -d "$sample_dir/mRNA" ]] && safe_ln "$sample_dir/mRNA" "$dest/mrna"

    # TCR/BCR folders (只链接文件夹，不再额外链接单个xlsx，避免重复入口)
    [[ -d "$sample_dir/TCR" ]] && safe_ln "$sample_dir/TCR" "$dest/tcr"
    [[ -d "$sample_dir/BCR" ]] && safe_ln "$sample_dir/BCR" "$dest/bcr"

  done < <(find "$summary" -mindepth 1 -maxdepth 1 -type d -print0)
done

say "Done. Example links:"
find "$DEST_BASE" -maxdepth 2 -type l | sed -n '1,60p'
