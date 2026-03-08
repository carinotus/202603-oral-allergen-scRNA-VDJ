#!/usr/bin/env bash
# ============================================================
# 00_bash_commands_all.sh
# 联川10x (mRNA + VDJ) 数据从压缩包 → 解压 → 统一软链接到 data/raw/10x/
# 适用于：本项目的 3W_PBS_OVA / 5W_PBS_HDM_OVA / 3W_PBS_HDM
#
# 使用方式：
#   1) 在项目根目录运行：bash scripts/00_bash_commands_all.sh
#   2) 或按需要分别运行子脚本
# ============================================================

set -euo pipefail
# -e: 任意命令失败立即退出（避免“错了还继续跑”）
# -u: 使用未定义变量就报错退出（避免空变量导致误删/误解压）
# -o pipefail: 管道中任一命令失败则整个管道失败

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "[INFO] Project root: $ROOT"

# 0) 目录准备
mkdir -p data/raw data/processed data/external
mkdir -p data/raw/10x
mkdir -p results/figures results/tables results/logs
mkdir -p metadata

echo "[STEP] 1) OVA/5W (LC-like) 解压 + 软链接"
bash scripts/00_bash_commands_lc_like.sh # 解压联川格式的测序数据并软链接到 data/raw/10x/<sample_id>/{mrna,tcr,bcr}

echo "[STEP] 2) 3W_HDM (split) 校验 + 解压 + 软链接"
bash scripts/00_bash_commands_hdm_split.sh # 校验3W_PBS_HDM样格式的的Summary.tar.gz，解压到 unpacked/，再软链接到 data/raw/10x/<sample_id>/{mrna,tcr,bcr,loom}

echo "[STEP] 3) 验收统一入口 data/raw/10x"
# tree 用于直观看目录层级（若服务器没装 tree，可 sudo apt-get install tree）
if command -v tree >/dev/null 2>&1; then
  tree -L 2 data/raw/10x | sed -n '1,200p'
else
  echo "[WARN] tree not found, using find instead:"
  find data/raw/10x -maxdepth 2 -type l | sed -n '1,80p'
fi

echo "[STEP] 4) 生成 metadata/samples_auto.csv（用于R批量读取）"
# 你已有这个脚本就直接调用
bash scripts/00_make_samples_auto.sh # 生成 metadata/samples_auto.csv
bash scripts/00_make_metadata_from_raw10x.sh # 生成 metadata/samples.csv（包含更多手动补充的字段）

echo "[DONE] All steps finished."