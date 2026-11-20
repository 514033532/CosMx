# Fig 6 细胞通讯

## 面板
a NC/DN 通讯强度矩阵 b 细胞对 LR 数 DE 热图 c 上下调 LR 对条形 d 中心度 DE 热图

## 上游脚本（scripts_archive/）
cellchat1-code.sh/cellchat1-stlearn-code.sh（LR 对箱线）、cellchat-cellpair/lrpair-热图code.sh + stlearn-cellpair-热图code.sh（三层热图）、cellchat1-stlearn-code-热图.sh；数据 comm_strength/lr_num_de/lr_value_de/centrality_de.csv。

## 一键运行
```bash
./run_fig6.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
