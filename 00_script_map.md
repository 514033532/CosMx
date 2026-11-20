# 脚本 → 大图 总映射表（整理日期 2026-10-22）

来源指令：zhiling-20261020.txt。27 个散落脚本（已归档至 `scripts_archive/`，重命名消除歧义）→ 每张大图一个目录、一个入口脚本。

## 图目录结构

```
code_organized_20261022/
├── 00_script_map.md          ← 本表
├── scripts_archive/          ← 27 个原始脚本（未改动）
├── lib/common.py             ← 共享样式/路径（OUT 指向本包）
├── source_data/ → 软链 NC_manuscript_v1_20261009/source_data
├── figures/                  ← 所有大图输出（PDF/PNG/TIFF）
├── fig1_atlas/ … fig7_mIF/ sfig1_NMF/ sfig2_QC/
│    ├── README.md            ← 本图面板↔脚本↔数据文件映射
│    ├── run_figN.sh          ← 一键入口（一个代码得全图）
│    └── plot_figN.py         ← 出图脚本（来自稿件包 code/）
```

## 27 脚本 → 大图映射总表

| # | 归档名 | 原路径 | 产出 | 归属面板 |
|---|--------|--------|------|----------|
| 1 | 0514.sh | YUBO-202508/0514.sh | zone 7/9/10 Signature top20 热图（DNvsNC + DN1DN2vsNC） | f2 区域分带 |
| 2 | 0516.sh | YUBO-202508/0516.sh | 同上，A4 论文版（zone 1/7/9/10） | f2（0514 的后续版） |
| 3 | CODE1-GROUP.sh | 顶层 | 细胞比例 DE 箱线（DN1 vs NC，逐细胞 PDF/PNG） | f4 比例箱线 |
| 4 | cellchat1-code.sh | 顶层 | cellchat LR 对通讯分数箱线（前 15 对） | f6 箱线 |
| 5 | CODE1-GROUP-NC-DN.sh | 顶层 | 细胞比例 DE 箱线（DN vs NC）+ 合并图 | f4 比例箱线 |
| 6 | CODE1-GROUP-指令.sh | 顶层 | **需求文档**（DN1 vs DN2 箱线规格） | 脚本 10 的需求来源 |
| 7 | NC-DN2-code-热图.sh | 顶层 | 共定位富集（DN2 vs NC）up/down top20 热图 | f4 共定位热图 |
| 8 | CODE1-GROUP_DN2_vs_NC.sh | 顶层 | 细胞比例 DE 箱线（DN2 vs NC） | f4 比例箱线 |
| 9 | cellchat1-stlearn-code.sh | 顶层 | stlearn LR 对箱线（前 15 对） | f6 箱线 |
| 10 | CODE1-GROUP_DN1_vs_DN2.sh | 顶层 | 细胞比例 DE 箱线（DN1 vs DN2） | f4 比例箱线 |
| 11 | 7个文件循环-code-热图.sh | 顶层 | 通用批量热图 v1（12 文件，20250826 数据） | f2/f4/f6 热图 |
| 12 | cellchat1-stlearn-code-热图.sh | 顶层 | stlearn lr_pair DE 热图 | f6 热图 |
| 13 | CODE1-GROUP-NC-DN-热图代码.sh | 顶层 | 细胞比例 DE 热图（文件名残留 LR_pair） | f4 比例热图 |
| 14 | 7个文件循环-code-热图-v2.sh | 顶层 | 批量热图 v2（22 文件，20251114 数据） | f2/f4/f6 热图 |
| 15 | 7个文件循环-code-热图-v2-BU.sh | 顶层 | v2 裁剪版：仅 centrality 3 类（新增矩阵文件处理） | f2 中心性热图 |
| 16 | 7个文件循环-code-热图-v2-BU2.sh | 顶层 | v2 裁剪版：仅 Moran + 细胞比例 | f2 moran / f4 比例 |
| 17 | 每个细胞类型的差异分析结果.sh | 顶层 | 11 细胞类型 Seurat DEG + GO + 热图（v1） | f5/f3 细胞 DEG |
| 18 | 每个细胞类型的差异分析结果2.sh | 顶层 | 同 17 的鲁棒修复版（top50、tryCatch、输出位置不同） | f5/f3（最终版） |
| 19 | ranse/code.sh | ranse/ | 一行 QuPath 路径入口 | f7 工具链 |
| 20 | run_all.sh | panel2_analysis/scripts/ | QuPath 批量：04_process_slide + 05_export_images | f7 切片导出 |
| 21 | CODE.deg.sh | 20251114_share/groups/cell_type/ | 8 细胞类型 edgeR 火山合并 4×2 | f5（或 f3）火山 |
| 22 | colocate-热图code.sh | 20250827/.../cell_colocate/ | 共定位 CellEnrich（DN2 vs NC）热图 | f4 共定位 |
| 23 | zone-code.sh | 20250827/.../cell_zone/ | zone 0/3/6 Signature t 检验 top20 热图 ×2 注释 | f2 zone 热图 |
| 24 | enrich-coe.sh | 20250827/.../Enrichment_results/ | zone 0/3/6 KEGG 气泡三联（v1） | f3 KEGG |
| 25 | enrich-kegg.sh | 同上 | 同 24 的最终版（名称换行、大字体）→ KEGG_Enrichment_Wrapped_Name | f3 KEGG（最终） |
| 26 | enrich-zone1-7-9-10.sh | 同上 | zone 1/7/9/10 KEGG 四联（输出未见成品，未运行） | f3 KEGG 补充 |
| 27 | cellchat-cellpair-热图code.sh | 20250827/.../cellchat/cell_pair/ | cellchat 细胞对 LR 数 DE top20 热图 ×2 注释 | f6 热图 |
| 28 | cellchat-lrpair-热图code.sh | 20250827/.../cellchat/lr_pair/ | cellchat LR 值 DE top20 热图 ×2 注释 | f6 热图 |
| 29 | stlearn-cellpair-热图code.sh | 20250827/.../stlearn/cell_pair/ | stlearn 细胞对 LR 数 DE 热图（与 27 同模板） | f6 热图 |

## 版本谱系

- 0514 → 0516（A4 版）；11(v1, 20250826) → 14(v2, 20251114) → 15(BU, centrality) / 16(BU2, moran+比例)
- 17 → 18（Seurat DEG 流程修复版）；4 ≈ 9（cellchat/stlearn 同模板）；3/5/8/10（同一箱线模板 × 4 组比较）
- 24 → 25（KEGG 气泡最终版）；27/28/29（同模板 × 3 数据源）
- f3 KEGG：zone 0/3/6（25）+ zone 1/7/9/10（26，未运行，需补跑或从 cell_zone/KEGG_Enrichment_Wrapped_Name 取）

## 已发现疑点（待用户确认）

1. **colocate-热图code.sh（#22）**：`setwd` 指向 Group1_DN2-VS-DN1，而读的是 Group1_DN2-VS-NC 数据——输出落位可能错乱。
2. **enrich-zone1-7-9-10.sh（#26）**：输出 `KEGG_Enrichment_CellZone` 未见成品。
3. 20250827 目录里 10-15 日与 10-20 日两批输出并存（新旧两版脚本），新图以 10-20 日的 `cell_communication_*_heatmap_group1/2` 为准。
4. 人工截图命名笔误：`1889-2.PNG`（疑 2889-2）、`8119.5.PNG`（疑 8119-5）——肾小球目录。
5. R 环境缺包：tidyr/ggpubr/stringr/Seurat/clusterProfiler/ggrepel——归档脚本如需重跑需先补装（conda 或 CRAN）。
