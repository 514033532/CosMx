# Fig 3 区域病理

## 面板
a-c zone 0/3/6 火山 + KEGG 通路条形

## 上游脚本（scripts_archive/）
enrich-coe.sh→enrich-kegg.sh（zone 0/3/6 KEGG 最终版）、enrich-zone1-7-9-10.sh（1/7/9/10，未运行）；每个细胞类型的差异分析结果2.sh（GO 部分）；数据 zone_deg/zone_pathways.csv。

## 一键运行
```bash
./run_fig3.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
