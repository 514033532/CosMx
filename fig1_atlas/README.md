# Fig 1 图谱总览

## 面板
a 流程图 b 大类UMAP c 18类UMAP d dotplot e 空间+16样本网格 f 比例堆叠

## 上游脚本（scripts_archive/）
无直接上游 R 脚本候选——数据来自 20251114_share/cell_anno/filter.cell_anno.h5ad（UMAP/细胞注释）与 sample_celltype_proportion.csv；经稿件包 extract_cosmx.py 出 umap.csv/cells.csv/dotplot_values.csv/major_group_proportions.csv。

## 一键运行
```bash
./run_fig1.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
