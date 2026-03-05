#!/usr/bin/env bash
# ============================================================
# 00_bash_commands_lc_like.sh
# 适用结构：
#   data/raw/3W_PBS_OVA/Summary.tar.gz
#   data/raw/5W_PBS_HDM_OVA/Summary.tar.gz
# 解压后通常为：
#   unpacked/Summary/<sample_id>/mRNA  (barcodes/features/matrix)
#   unpacked/Summary/<sample_id>/TCR, BCR (filtered_contig_annotations.xlsx 等)
#
# 作用：
#   1) 找压缩包 / 预览结构
#   2) 解压到 unpacked/
#   3) 用已有脚本 00_link_10x_lc.sh 把样本软链接到 data/raw/10x/<sample_id>/{mrna,tcr,bcr}
# ============================================================

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# ---- 项目里的实际位置 ----
D1="data/raw/3W_PBS_OVA"
D2="data/raw/5W_PBS_HDM_OVA"

echo "[lc_like] 1) 查找压缩包（tar.gz/tgz/zip）"
# find 参数解释：
# -maxdepth 2：只搜到第二层，避免把整个目录树扫爆
# -type f：只找文件
# \( ... \)：把多个 -name 条件用 OR 连起来
find "$D1" "$D2" -maxdepth 2 -type f \
  \( -name "*.tar.gz" -o -name "*.tgz" -o -name "*.tar" -o -name "*.zip" \) -print || true

# 假设联川提供的是 Summary.tar.gz（如果文件名不同，改这里）
F1="$D1/Summary.tar.gz"
F2="$D2/Summary.tar.gz"

echo "[lc_like] 2) 预览压缩包内部结构（不解压）"
# tar -t：list 内容；-f 指定文件；head -n 40 只看前40行
if [[ -f "$F1" ]]; then
  tar -tf "$F1" | head -n 40
else
  echo "[WARN] not found: $F1"
fi
if [[ -f "$F2" ]]; then
  tar -tf "$F2" | head -n 40
else
  echo "[WARN] not found: $F2"
fi

echo "[lc_like] 3) 创建解压目录（-p 表示父目录不存在也一并创建，不报错）"
mkdir -p "$D1/unpacked" "$D2/unpacked"

echo "[lc_like] 4) 解压（tar -x 解包；-z 处理gz；-f 指定文件；-C 指定解压目标目录）"
if [[ -f "$F1" ]]; then
  tar -xzf "$F1" -C "$D1/unpacked"
fi
if [[ -f "$F2" ]]; then
  tar -xzf "$F2" -C "$D2/unpacked"
fi

echo "[lc_like] 5) 找 mRNA / TCR / BCR（用于确认解压正确）"
# -type d 找目录；-name "mRNA" 按名字匹配
find "$D1/unpacked" "$D2/unpacked" -maxdepth 6 -type d -name "mRNA" -print | sed -n '1,50p'
find "$D1/unpacked" "$D2/unpacked" -maxdepth 8 -type d -name "TCR" -print | sed -n '1,50p'
find "$D1/unpacked" "$D2/unpacked" -maxdepth 8 -type d -name "BCR" -print | sed -n '1,50p'

echo "[lc_like] 6) 统一软链接到 data/raw/10x（调用你已有脚本）"
# 这个脚本会把 unpacked/Summary/<sample>/mRNA,TCR,BCR 链接到 data/raw/10x/<sample>/
bash scripts/00_link_10x_lc.sh

echo "[lc_like] done."