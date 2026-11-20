# Fig 5 三维差异

## 面板
a 比例维度 b 聚类维度 c 中心度维度 火山

## 上游脚本（scripts_archive/）
CODE.deg.sh（8 细胞类型火山合并）、每个细胞类型的差异分析结果2.sh（11 细胞类型 DEG+GO 热图）；数据 ct_deg/centrality_de.csv。

## 一键运行
```bash
./run_fig5.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
