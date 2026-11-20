# SFig 1 NMF

## 面板
a MP 热图 b MP 箱线 c 通路气泡

## 上游脚本（scripts_archive/）
数据 gene_NMF/cNMF.metadata.csv、mean_MPs_sample.csv、MP_C2_pathway.csv（20251114_share）。

## 一键运行
```bash
./run_sfig1.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
