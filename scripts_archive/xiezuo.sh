fig1 (zhubiaoti)
(A)CosMx SMI基因技术原理

b,样本细胞空间原位展示he cell  type zhushi


c,细胞大类UMAP 根据细胞聚类分析结果, 对每个Cluster进行细胞大类注释; 再分别对不同细胞大类单独进行聚类并注释细胞亚类。

d,细胞亚类UMAP


e结果说明: 该模块展示整合所有样本后的细胞UMAP图谱, 并对样本层面(Sample)以及聚类层面(Cluster)的信息进行Overview展示。

f细胞大类Marker基因分析
g各样本中细胞类型的分布




fig2 (zhubiaoti)

a,空间结构域鉴定结果原位展示;
b,空间结构域UMAP;
c,zone_cell_proportion展示了每个结构域中细胞类型的分布。对于空间结构域的注释并没有统一的标准, 一般而言,如果空间结构域与病理组织结构相对吻合, 可以用病理组织结构信息来注释不同空间结构域;
d,每个结构域的Marker基因;
e,sample_zone_proportion展示了每个样本中不同结构域的比例分布;
f,_DN1_DN2-VS-NC空间结构域的比例组间差异分析;
g,每个样本不同空间结构域之间的关联性差异分析;




fig3 (zhubiaoti)
a,0 3 6 cell_zone/Enrichment_results/不同空间结构域的组间差异基因富集分析
b,0 3 6 cell_zone不同空间结构域的组间差异基因热图
c,0 3 6 cell_zone 不同空间结构域的组间差异基因火山图



fig4 (zhubiaoti)
a,spatial_distr/centrality，cell_closeness_heatmap_group1.png细胞层面空间分布特征的组间差异分析
b,不同组间的基因空间分布特征差异与组间差异基因不同, 前者反应的是空间特征的分布变化, 后者只是表达水平的变化。
c,不同细胞类型邻近空间富集特征差异分析
d,基于CellChat的组间细胞互作差异cell_pair
e,基于CellChat的组间细胞互作差异lr_pair，
f,基于stlearn的组间细胞互作差异cell_pair
g,基于stlearn基于CellChat的组间细胞互作差异lr_pair，


supplementary fig1 (zhubiaoti)
a,基于cNMF分析获得的基因表达程序（MP）

b,组间差异MP分析
c,每个MP的KEGG


fig5 (zhubiaoti)
### 图1各子图标题提取（A、B、C、D、F）
- **A**：CellCloseness_DE_list Expression (DN1, DN2 vs NC)
- **B**：gene_moran Expression (DN1, DN2 vs NC)
- **C**：cell_enrich Expression (DN1, DN2 vs NC)
- **D**：cell_communication cellchat cellpair Expression (DN1, DN2 vs NC)
 **e**：cell_communication_cellchat_cell lrpair Expression (DN1, DN2 vs NC)
- **F**：cell_communication_stlearn_cellpair Expression (DN1, DN2 vs NC)

 **g**：cell_communication_stlearn_lrpair Expression (DN1, DN2 vs NC)


fig6 (zhubiaoti)

### 图1各子图标题提取（整合两张图的子图信息）
- **A**：CellType_prop_DE_list Expression (DN1, DN2 vs NC)
- **B**：CellClustering_DE_list Expression (DN1, DN2 vs NC)
- **C**：CellDegree_DE_list Expression (DN1, DN2 vs NC)

e,/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/Cortical_fibroblasts

f,/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/M2_Macrophage

g,/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/Monocyte_cell

h,/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/T
i,/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/B








### 各热图差异细胞类型（up/down in case）提取

#### 1. **CellType_prop_DE_list Expression (DN1, DN2 vs NC)**
- **Up in case（红色/粉色）**：
  Mast, T, B, NK, Plasma, Cortical.fibroblasts, M2.Macrophage, Monocyte.cell
- **Down in case（绿色）**：
  Podocytes, Proximal.tubular.cell, Distal.tubular.cell

---

#### 2. **CellClustering_DE_list Expression (DN1, DN2 vs NC)**
- **Up in case（红色/粉色）**：
  Plasma, T, Activated.myofibroblasts
- **Down in case（绿色）**：
  Monocyte.cell, M2.Macrophage, Cortical.fibroblasts, Distal.tubular.cell

---

#### 3. **CellDegree_DE_list Expression (DN1, DN2 vs NC)**
- **Up in case（红色/粉色）**：
  Cortical.fibroblasts, M2.Macrophage, Plasma, Activated.myofibroblasts
- **Down in case（绿色）**：
  Monocyte.cell, Proximal.tubular.cell, Distal.tubular.cell

---

如果你需要，我可以把这些结果整理成一个**对比表格**，方便你直接放到论文里。需要吗？

congshangdaoxia,zuodaoyou,
/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/Cortical_fibroblasts

/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/M2_Macrophage

/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/Monocyte_cell

/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/T
/home/zhh2/YUBO-202508/20251114_share/groups/cell_type/DEG_result/Heatmap/B








### 一、中文版本（符合论文发表规范）

---

#### **Fig 1. 糖尿病肾病（DKD）的单细胞及空间转录组细胞图谱构建与注释**
- **(A)** CosMx SMI空间转录组技术原理示意图
- **(B)** 样本细胞空间原位分布及细胞类型注释
- **(C)** 基于无监督聚类细胞大类UMAP聚类及注释
- **(D)** 细胞亚类UMAP聚类及注释
- **(E)** 整合所有样本的细胞UMAP图谱
- **(F)** 细胞大类Marker基因表达特征分析
- **(G)** 各样本中细胞类型的分布比例

---

#### **Fig 2. 糖尿病肾病肾脏空间结构域的鉴定与特征分析**
- **(A)** 空间结构域鉴定结果的原位可视化
- **(B)** 空间结构域UMAP聚类分析
- **(C)** 各空间结构域中细胞类型的分布比例
- **(D)** 各空间结构域的Marker基因表达分析
- **(E)** 各样本中不同空间结构域的比例分布
- **(F)** DN1/DN2 vs NC空间结构域比例的组间差异分析
- **(G)** 各样本内不同空间结构域之间的关联性差异分析

---

#### **Fig 3. 空间结构域的组间差异基因及功能富集分析（DN1, DN2 vs NC）**
- **(A)** 不同空间结构域的组间差异基因功能富集分析
- **(B)** 不同空间结构域的组间差异基因表达热图
- **(C)** 不同空间结构域的组间差异基因火山图

---

#### **Fig 4. 细胞空间分布与细胞互作的组间差异特征分析**
- **(A)** 细胞层面空间邻近度（Cell closeness）的组间差异分析
- **(B)** 基因空间自相关性（Moran’s I）的组间差异分析
- **(C)** 不同细胞类型的邻近空间富集特征差异分析
- **(D)** 基于CellChat的组间细胞互作差异（Cell pair水平）
- **(E)** 基于CellChat的组间细胞互作差异（配体-受体对水平）
- **(F)** 基于stlearn的组间细胞互作差异（Cell pair水平）
- **(G)** 基于stlearn的组间细胞互作差异（配体-受体对水平）

---

#### **Supplementary Fig 1. 基因表达程序（Metaprograms）的鉴定与功能分析**
- **(A)** 基于cNMF分析鉴定的基因表达程序（Metaprograms, MPs）
- **(B)** 组间差异基因表达程序（MPs）分析
- **(C)** 各基因表达程序（MPs）的KEGG通路富集分析

---

#### **Fig 5. 空间与互作特征的差异表达分析（DN1, DN2 vs NC）**
- **(A)** 细胞邻近度差异表达分析（DN1, DN2 vs NC）
- **(B)** 基因空间自相关性差异表达分析（DN1, DN2 vs NC）
- **(C)** 细胞邻近空间富集差异表达分析（DN1, DN2 vs NC）
- **(D)** 基于CellChat的细胞互作差异表达分析（Cell pair水平，DN1, DN2 vs NC）
- **(E)** 基于CellChat的细胞互作差异表达分析（配体-受体对水平，DN1, DN2 vs NC）
- **(F)** 基于stlearn的细胞互作差异表达分析（Cell pair水平，DN1, DN2 vs NC）
- **(G)** 基于stlearn的细胞互作差异表达分析（配体-受体对水平，DN1, DN2 vs NC）

---

#### **Fig 6. 关键细胞亚群的组间差异表达分析（DN1, DN2 vs NC）**
- **(A)** 细胞类型比例差异表达分析（DN1, DN2 vs NC）
- **(B)** 细胞聚类差异表达分析（DN1, DN2 vs NC）
- **(C)** 细胞中心度差异表达分析（DN1, DN2 vs NC）
- **(D)** 皮质成纤维细胞（Cortical fibroblasts）的组间差异基因热图
- **(E)** M2型巨噬细胞（M2 Macrophage）的组间差异基因热图
- **(F)** 单核细胞（Monocyte cell）的组间差异基因热图
- **(G)** T细胞（T cells）的组间差异基因热图
- **(H)** B细胞（B cells）的组间差异基因热图

---

### 二、英文版本（符合国际期刊发表规范）

---

#### **Fig 1. Construction and Annotation of Single-Cell and Spatial Transcriptomic Cell Atlas in Diabetic Kidney Disease (DKD)**
- **(A)** Schematic diagram of the CosMx SMI spatial transcriptomics technology
- **(B)** Spatial in situ visualization of cells and cell type annotation in samples
- **(C)**  Based on unsupervised clustering UMAP clustering and annotation of major cell types
- **(D)** UMAP clustering and annotation of cell subtypes
- **(E)** Integrated UMAP map of cells across all samples
- **(F)** Marker gene expression analysis of major cell types
- **(G)** Distribution proportion of cell types across samples

---

#### **Fig 2. Identification and Characterization of Spatial Domains in Diabetic Kidney Disease**
- **(A)** In situ visualization of identified spatial domains
- **(B)** UMAP clustering analysis of spatial domains
- **(C)** Distribution proportion of cell types within each spatial domain (spatial domains are annotated based on pathological tissue structures)
- **(D)** Marker gene expression analysis of each spatial domain
- **(E)** Proportional distribution of different spatial domains across samples
- **(F)** Intergroup differential analysis of spatial domain proportions in DN1/DN2 vs NC
- **(G)** Correlation difference analysis between different spatial domains within each sample

---

#### **Fig 3. Intergroup Differentially Expressed Genes and Functional Enrichment Analysis of Spatial Domains (DN1, DN2 vs NC)**
- **(A)** Functional enrichment analysis of intergroup differentially expressed genes (DEGs) across spatial domains
- **(B)** Expression heatmap of intergroup DEGs across spatial domains
- **(C)** Volcano plot of intergroup DEGs across spatial domains

---

#### **Fig 4. Intergroup Differential Analysis of Cell Spatial Distribution and Cell-Cell Communication**
- **(A)** Intergroup differential analysis of cell-level spatial closeness
- **(B)** Intergroup differential analysis of gene spatial autocorrelation (Moran’s I)
- **(C)** Differential analysis of adjacent spatial enrichment characteristics across cell types
- **(D)** Intergroup cell-cell communication differences based on CellChat (Cell pair level)
- **(E)** Intergroup cell-cell communication differences based on CellChat (Ligand-receptor pair level)
- **(F)** Intergroup cell-cell communication differences based on stlearn (Cell pair level)
- **(G)** Intergroup cell-cell communication differences based on stlearn (Ligand-receptor pair level)

---



### 英文极简发表版（最短、清晰、符合期刊格式）

#### **Fig 2. Identification and characterization of spatial domains in diabetic kidney disease**
- (A) In situ visualization of spatial domains
- (B) UMAP clustering of spatial domains
- (C) Cell type proportion in each spatial domain
- (D) Marker gene expression of spatial domains
- (E) Spatial domain distribution across samples
- (F) Proportional differences of spatial domains (DN1/DN2 vs NC)
- (G) Correlation between spatial domains within samples

---

#### **Fig 3. Differentially expressed genes and functional enrichment of spatial domains (DN1, DN2 vs NC)**
- (A) Functional enrichment of DEGs in spatial domains
- (B) Heatmap of DEGs in spatial domains
- (C) Volcano plot of DEGs in spatial domains

---

#### **Fig 4. Differences in cell spatial distribution and cell–cell communication**
- (A) Cell spatial closeness analysis
- (B) Gene spatial autocorrelation (Moran’s I)
- (C) Spatial enrichment of cell types
- (D) Cell–cell communication by CellChat (cell pair)
- (E) Ligand–receptor pairs by CellChat
- (F) Cell–cell communication by stlearn (cell pair)
- (G) Ligand–receptor pairs by stlearn












#### **Supplementary Fig 1. Identification and Functional Analysis of Gene Expression Programs (Metaprograms)**
- **(A)** Gene expression programs (Metaprograms, MPs) identified by cNMF analysis
- **(B)** Intergroup differential analysis of gene expression programs (MPs)
- **(C)** KEGG pathway enrichment analysis of each gene expression program (MP)

---

#### **Fig 5. Differential Expression Analysis of Spatial and Communication Features (DN1, DN2 vs NC)**
- **(A)** Differential expression analysis of cell closeness (DN1, DN2 vs NC)
- **(B)** Differential expression analysis of gene spatial autocorrelation (DN1, DN2 vs NC)
- **(C)** Differential expression analysis of cell adjacent spatial enrichment (DN1, DN2 vs NC)
- **(D)** Expression analysis of intergroup cell-cell communication differences based on CellChat (Cell pair level, DN1, DN2 vs NC)
- **(E)** Expression analysis of intergroup cell-cell communication differences based on CellChat (Ligand-receptor pair level, DN1, DN2 vs NC)
- **(F)** Expression analysis of intergroup cell-cell communication differences based on stlearn (Cell pair level, DN1, DN2 vs NC)
- **(G)** Expression analysis of intergroup cell-cell communication differences based on stlearn (Ligand-receptor pair level, DN1, DN2 vs NC)

---### 中文（极简·发表版）
#### Fig 5. 空间特征与细胞互作差异表达分析（DN1、DN2 vs NC）
- (A) 细胞邻近度差异表达分析
- (B) 基因空间自相关差异表达分析
- (C) 细胞空间邻近富集差异分析
- (D) 基于CellChat的细胞对互作差异分析
- (E) 基于CellChat的配体-受体对互作差异分析
- (F) 基于stlearn的细胞对互作差异分析
- (G) 基于stlearn的配体-受体对互作差异分析

---

### English (Concise & Publication-ready)
#### Fig 5. Differential Expression Analysis of Spatial and Cell-Cell Communication Features (DN1, DN2 vs NC)
- (A) Differential expression analysis of cell closeness
- (B) Differential expression analysis of gene spatial autocorrelation
- (C) Differential analysis of cell spatial adjacent enrichment
- (D) Cell–cell communication analysis by CellChat (cell pair)
- (E) Ligand–receptor pair analysis by CellChat
- (F) Cell–cell communication analysis by stlearn (cell pair)
- (G) Ligand–receptor pair analysis by stlearn

#### **Fig 6. Intergroup Differential Expression Analysis of Key Cell Subpopulations (DN1, DN2 vs NC)**
- **(A)** Differential expression analysis of cell type proportions (DN1, DN2 vs NC)
- **(B)** Differential expression analysis of cell clustering (DN1, DN2 vs NC)
- **(C)** Differential expression analysis of cell centrality (DN1, DN2 vs NC)
- **(D)** Heatmap of intergroup DEGs in Cortical fibroblasts
- **(E)** Heatmap of intergroup DEGs in M2 Macrophage
- **(F)** Heatmap of intergroup DEGs in Monocyte cell
- **(G)** Heatmap of intergroup DEGs in T cells
- **(H)** Heatmap of intergroup DEGs in B cells

---#### **Fig 6. Differential expression analysis of key cell subsets (DN1, DN2 vs NC)**
- (A) Cell type proportion analysis
- (B) Cell clustering analysis
- (C) Cell centrality analysis
- (D) DEG heatmap of cortical fibroblasts
- (E) DEG heatmap of M2 macrophages
- (F) DEG heatmap of monocytes
- (G) DEG heatmap of T cells
- (H) DEG heatmap of B cells

如果你需要，我可以进一步帮你把每个子图的**图例（legend）和统计方法说明**也补充成中英文版本，直接适配论文格式。需要我补充吗？





