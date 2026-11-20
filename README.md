# CosMx SMI analysis code — human diabetic kidney disease

Analysis and figure-generation code for:

**Single-cell-resolution spatial profiling of human diabetic kidney disease reveals remodelled intercellular communication, immune-enriched domains and reduced glomerular CTLA4**

Zongchao Yu, Huanhuan Lin, Liangyan Wu, Wenyu Gong, Lingyu Huang, Xiaoni Zhong, Meiqian Li, Zhipeng Zeng, Yangna Shi, Han Zhang, Jingsheng Ma, Donge Tang

## Repository structure

```
code_organized_20261022/
├── 00_script_map.md          # mapping of archived scripts → manuscript figures
├── PACKAGE_VERSIONS.md       # environment, package versions, random seeds
├── README.md
├── lib/common.py             # shared styles/paths
├── source_data/              # per-panel source CSVs (see Data below)
├── figures/                  # generated figure files (PDF/PNG/TIFF)
├── fig1_atlas/ … fig7_mIF/   # one directory per figure
│   ├── README.md             # panel ↔ script ↔ data mapping for that figure
│   ├── run_figN.sh           # one-command entry point
│   └── plot_figN.py          # plotting script
├── sfig1_NMF/ sfig2_QC/      # supplementary figures
├── fig8.py                   # additional panel
└── scripts_archive/          # 27 original scripts (historical, unmodified)
```

## Requirements

- **R**: conda environment `R4.5` — R 4.5.1 with Seurat 5.5.1, CellChat 1.6.1, edgeR 4.8.2, fgsea 1.36.2, clusterProfiler 4.18.4, ClusterGVis 0.1.4 and tidyverse/plotting packages.
- **Python**: 3.13 with scanpy 1.12.4, pandas 2.3.3, numpy 2.4.1, scipy 1.18.1, matplotlib 3.10.8, anndata 0.13.4, Pillow.
- Full version table: [`PACKAGE_VERSIONS.md`](PACKAGE_VERSIONS.md).

## Quick start

Each figure directory has a one-command entry point, e.g.:

```bash
cd fig1_atlas && ./run_fig1.sh      # outputs ../figures/Fig1_*.{pdf,png,tif}
```

Figure outputs are written to `figures/`.

## Data

- **`source_data/`**: per-panel source CSV files used by the plotting scripts (sample- and cell-level summaries). The largest file is ~89 MB.
- **Processed CosMx data**: deposited in the OMIX database (CNCB-NGDC), accession **OMIX018638** (https://ngdc.cncb.ac.cn/omix/preview/GzaeetoY7K).
- **mIF crop images** (Fig 7 panel B): provided outside this repository. Set the environment variable `MIF_CROPS_DIR` to the folder containing per-slide crop subdirectories before running `fig7_mIF`; by default the scripts look next to the analysis workspace.
- **Supplementary tables**: `Supplementary_Tables.xlsx` in the manuscript package.

## Reproducibility notes

- Random seeds used in plotting scripts are listed in [`PACKAGE_VERSIONS.md`](PACKAGE_VERSIONS.md).
- Metaprogram (cNMF) and stLearn summaries were computed upstream from the raw CosMx data; the upstream environment is documented in [`PACKAGE_VERSIONS.md`](PACKAGE_VERSIONS.md).

## License

MIT — see [LICENSE](LICENSE).
