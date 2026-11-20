# Fig 4 分组重塑

## 面板
a 比例DE热图 b 基因模块热图 c 共定位DE热图 d 代表对箱线

## 上游脚本（scripts_archive/）
CODE1-GROUP.sh/CODE1-GROUP-NC-DN.sh/CODE1-GROUP_DN2_vs_NC.sh/CODE1-GROUP_DN1_vs_DN2.sh（4 组比较箱线）、CODE1-GROUP-NC-DN-热图代码.sh+BU2（比例热图）、NC-DN2-code-热图.sh/colocate-热图code.sh（共定位热图）；数据 ct_prop_de/enrich_de/enrich_zscores/gene_module_*.csv。

## 一键运行
```bash
./run_fig4.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
