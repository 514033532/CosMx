# SFig 2 QC 稳健性

## 面板
a QC 概览 b 阈值敏感性 c 共表达重叠

## 上游脚本（scripts_archive/）
数据 qc_overview/mif_sensitivity/mif_coexpression.csv；QuPath qc 日志（run_all.sh）。

## 一键运行
```bash
./run_sfig2.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
