#!/usr/bin/env bash
# ============================================================
# 00_bash_commands_hdm_split.sh
# 适用结构（你截图那种）：
#   data/raw/3W_PBS_HDM/mrnaseq/Summary.tar.gz (+ Summary.tar.gz.md5)
#   data/raw/3W_PBS_HDM/tcrseq/Summary.tar.gz  (+ Summary.tar.gz.md5)
#   data/raw/3W_PBS_HDM/bcrseq/Summary.tar.gz  (+ Summary.tar.gz.md5)
#
# 解压后通常为：
#   mrnaseq/unpacked/Summary/1_Cellranger_result/<sample>/filtered_feature_bc_matrix/
#   tcrseq/unpacked/Summary/1_Cellranger_result/<sample>/...
#   bcrseq/unpacked/Summary/1_Cellranger_result/<sample>/...
#
# 作用：
#   1) md5 校验（可选但强烈建议）
#   2) 分别解压到各自 unpacked/
#   3) 调用 00_link_10x_3W_PBS_HDM.sh 软链接到 data/raw/10x/<sample>/{mrna,tcr,bcr,loom}
# ============================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BASE="data/raw/3W_PBS_HDM"
MRNA_DIR="$BASE/mrnaseq"
TCR_DIR="$BASE/tcrseq"
BCR_DIR="$BASE/bcrseq"

echo "[hdm_split] 1) （可选）md5 校验：md5sum -c 会读取 .md5 文件逐行校验"
# 参数解释：
# md5sum -c xxx.md5：按文件里记录的md5逐个检查，输出 OK/FAILED
# 用 (cd dir && cmd) 避免写一堆长路径
if [[ -f "$MRNA_DIR/Summary.tar.gz.md5" ]]; then (cd "$MRNA_DIR" && md5sum -c Summary.tar.gz.md5); fi
if [[ -f "$TCR_DIR/Summary.tar.gz.md5"  ]]; then (cd "$TCR_DIR"  && md5sum -c Summary.tar.gz.md5); fi
if [[ -f "$BCR_DIR/Summary.tar.gz.md5"  ]]; then (cd "$BCR_DIR"  && md5sum -c Summary.tar.gz.md5); fi

echo "[hdm_split] 2) 创建解压目录"
mkdir -p "$MRNA_DIR/unpacked" "$TCR_DIR/unpacked" "$BCR_DIR/unpacked"

echo "[hdm_split] 3) 解压三个 Summary.tar.gz 到各自 unpacked/"
# tar -xzf：解压 tar.gz；-C：指定输出目录
if [[ -f "$MRNA_DIR/Summary.tar.gz" ]]; then tar -xzf "$MRNA_DIR/Summary.tar.gz" -C "$MRNA_DIR/unpacked"; fi
if [[ -f "$TCR_DIR/Summary.tar.gz"  ]]; then tar -xzf "$TCR_DIR/Summary.tar.gz"  -C "$TCR_DIR/unpacked";  fi
if [[ -f "$BCR_DIR/Summary.tar.gz"  ]]; then tar -xzf "$BCR_DIR/Summary.tar.gz"  -C "$BCR_DIR/unpacked";  fi

echo "[hdm_split] 4) 解压后快速验证：是否存在 filtered_feature_bc_matrix"
find "$MRNA_DIR/unpacked" -maxdepth 8 -type d -name "filtered_feature_bc_matrix" -print | sed -n '1,20p'

echo "[hdm_split] 5) 统一软链接到 data/raw/10x（调用你已有脚本）"
bash scripts/00_link_10x_3W_PBS_HDM.sh

echo "[hdm_split] done."