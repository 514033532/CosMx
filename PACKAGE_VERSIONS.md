# 软件包版本记录

**记录环境**：本工作站 miniconda（分析机），采集日期 2026-11-20。
**R 环境**：conda env `R4.5`（R version 4.5.1, 2025-06-13）
**Python 环境**：conda `base`（Python 3.13.2）

## R（conda env: R4.5）

| 包 | 版本 |
|---|---|
| Seurat | 5.5.1 |
| edgeR | 4.8.2 |
| fgsea | 1.36.2 |
| CellChat | 1.6.1（与正文一致 ✓） |
| GeneNMF | 未安装 —— 实际使用的元程序工具为 cNMF（Python，上游预计算），见下方备注 |
| clusterProfiler | 4.18.4 |
| org.Hs.eg.db | 3.22.0 |
| ClusterGVis | 0.1.4 |
| ggplot2 | 4.0.3 |
| ggpubr | 1.0.0 |
| ggsci | 5.2.0 |
| ggrepel | 0.9.8 |
| pheatmap | 1.0.13 |
| RColorBrewer | 1.1.3 |
| patchwork | 1.3.2 |
| cowplot | 1.2.0 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| tibble | 3.3.1 |
| stringr | 1.6.0 |
| reshape2 | 1.4.5 |
| stats | 4.5.1 |

## Python（conda base，Python 3.13.2）

| 包 | 版本 |
|---|---|
| python | 3.13.2 |
| scanpy | 1.12.4 |
| pandas | 2.3.3 |
| numpy | 2.4.1 |
| scipy | 1.18.1 |
| matplotlib | 3.10.8 |
| anndata | 0.13.4 |
| squidpy | 本机未安装（正文已不再声称使用） |
| stlearn | 本机未安装 —— Fig 5F/G 的 stLearn 通讯汇总为上游预计算，运行环境版本待补 |
| cnmf | 本机未安装 —— MP1–MP10 来自上游预计算数据 `20251114_share/gene_NMF/cNMF.metadata.csv`，运行环境版本待补 |

## 随机种子

自动扫描 `code_organized_20261022/` 的结果（2026-11-20）：

| 脚本 | 种子 |
|---|---|
| `fig8.py` | `np.random.seed(42)`；另一处 `np.random.default_rng(7)` |
| `fig4_group_remodelling/plot_fig4.py` | `np.random.seed(42)` |
| `sfig1_NMF/plot_sfig1.py` | `np.random.default_rng(42)` |
| `fig7_mIF/fig7_panels/panel_C.py` | `np.random.seed(42)` |
| `fig6_communication/fig6_panels/panel_E.py` | `np.random.default_rng(42)` |

以上多为出图抖动（jitter）种子。**Harmony / Leiden / cNMF 等分析步骤**若在原始分析脚本（可能不在本包内）中未显式设置种子，请在此写明，例如：

- Harmony/Leiden/cNMF 等步骤：未显式设置随机种子（使用各包默认值）
- （若设置了，请补脚本路径与数值）

## 待补（投稿前完成）

1. stLearn 与 cNMF 的上游运行环境及版本（若原始数据来自服务商/其他服务器，注明即可）；
2. Harmony / Leiden 的随机种子设置情况（见上）。
