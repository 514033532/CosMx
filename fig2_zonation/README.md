# Fig 2 功能分带

## 面板
a 代表样本 zone 空间图 b zone 组成 c zone marker dotplot d zone 比例 DE 热图 e 每 zone DEG 数

## 上游脚本（scripts_archive/）
scripts_archive/0514.sh→0516.sh（zone 1/7/9/10 Signature 热图）、zone-code.sh（zone 0/3/6 top20 热图）、7个文件循环-code-热图-v2-BU.sh（centrality）；数据经 extract 出 zone_markers/zone_composition/zone_prop_de/zone_deg.csv。

## 一键运行
```bash
./run_fig2.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
