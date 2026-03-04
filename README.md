# 202603 围断奶期口服过敏原免疫细胞图谱（10x mRNA + VDJ）

## 项目目标
- 系统解析围断奶期口服过敏原（HDM/OVA）在 3 周龄与 5 周龄条件下的免疫细胞图谱变化
- 联合 mRNA 状态与 VDJ 克隆扩增/谱系信息，解释关键细胞群的扩增与功能重塑

## 数据批次
- 3w HDM（mRNA + VDJ）
- 3w OVA（mRNA + VDJ）
- 5w HDM + OVA（mRNA + VDJ）

## 目录约定
- data/raw/：原始数据（只读，尽量软链接/不改动）
- data/processed/：中间对象（Seurat/VDJ 里程碑对象）
- metadata/：样本表 samples.csv、分组信息等
- scripts/：按步骤编号的可重复脚本
- R/：通用函数（读入、QC、作图、导出）
- results/figures/：图片输出（不入 git）
- results/tables/：最终表格输出（可选择入 git）
- results/logs/：每次运行日志（不入 git）

## 分析流程
01_import_qc.R
- 输入：10x mRNA 矩阵（h5 或 filtered_feature_bc_matrix）
- 输出：每样本 QC 图（nCount/nFeature/percent.mt 等）、过滤后的 Seurat 对象（里程碑1）

02_soup_decontam.R
- 输入：里程碑1对象 + 原始矩阵（raw/filtered）
- 输出：去环境RNA后对象（里程碑2）+ soup 估计报告

03_doublet.R
- 输入：里程碑2对象
- 输出：双细胞标注/去除后对象（里程碑3）+ 双细胞统计表

04_integration.R
- 输入：所有样本里程碑3对象
- 输出：整合对象（里程碑4）+ 批次效应评估图

05_annotation.R
- 输入：里程碑4
- 输出：注释后的对象（里程碑5）+ marker 表

06_vdj_analysis.R
- 输入：VDJ contig/clonotype 文件 + 里程碑5对象
- 输出：克隆扩增/多样性/共享克隆统计 + 图表

99_figures_export.R
- 输入：各里程碑对象与统计表
- 输出：最终图与表（results/figures, results/tables）
