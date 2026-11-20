# Fig 7 mIF 验证

## 面板
a 全景复合 b 代表肾小球 c 密度箱线 d 全片倍数 lollipop

## 上游脚本（scripts_archive/）
ranse/code.sh + run_all.sh（QuPath 04_process_slide/05_export_images）；人工截图 result20261020/{肾小球,肾小管间质,perivascular} 62 张（任务3将用于重排）；数据 mif_*.csv。

## 一键运行
```bash
./run_fig7.sh     # 输出: ../figures/<图名>.pdf/png/tif
```
依赖：python3 + matplotlib/pandas/PIL（同稿件包环境）；上游 R 片段重跑需 R 4.5.2 + pheatmap/ggplot2/dplyr/tidyr/ggpubr/stringr/tibble/RColorBrewer（缺 4 包待补装）。
