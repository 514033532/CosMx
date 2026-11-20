"""Extract CosMx source_data from 20251114_share (read-only). Run once; idempotent.
Schema of raw files probed 2026-10-09 (see ledger):
- DE list CSVs: ["", W, Pvalue, Significance, NC_Mean, <case>_Mean, Mean_difference, Change, 16 sample cols]
  entity names use dots ("Cortical.fibroblasts"), zones "X0".."X10", case-mean col varies per comparison.
- edgeR xlsx: [GeneSymbol, logFC, logCPM, PValue, FDR, Change] (files: DE_list_GroupX中*.xlsx)
- KEGG enrichment xlsx: [ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge, core_enrichment]
  (files: GroupX中*_GSEA_KEGG_enrichment.xlsx per zone dir)
- h5ad obs: sample/main_celltype(8 groups)/new_celltype(18+Mix)/zone_leiden; obsm banksy_umap; X float64 dense (chunked reads).
- cNMF.metadata.csv: 49,291 rows; main_celltype(8)/new_celltype(18)/zone_leiden(0-10 int)/MP1-3,5-10; x_slide_mm/y_slide_mm.
"""
import glob as gb
from pathlib import Path
import numpy as np
import pandas as pd
import common
from common import SD, SRC

COMP_DIRS = ["Group1_DN1-VS-NC", "Group1_DN2-VS-DN1", "Group1_DN2-VS-NC", "Group2_DN-VS-NC"]
GROUPS_MAP = {"Glomerular Cell": "Glomerular", "Interstitial Cell": "Interstitial",
              "Lymphoid cell": "Lymphoid", "Macrophage": "Macrophage",
              "Myeloid cell": "Myeloid", "Plasma": "Plasma",
              "Tubular Epithelial Cell": "Tubular", "Vascular Endothelial Cell": "Vascular"}
XLSX_MAP = {"Group1_DN1-VS-NC": "Group1中DN1比NC的差异",
            "Group1_DN2-VS-DN1": "Group1中DN2比DN1的差异",
            "Group1_DN2-VS-NC": "Group1中DN2比NC的差异",
            "Group2_DN-VS-NC": "Group2中DN比NC的差异"}
MEAN_COLS = {"Group1_DN1-VS-NC": ("NC_Mean", "DN1_Mean"),
             "Group1_DN2-VS-NC": ("NC_Mean", "DN2_Mean"),
             "Group1_DN2-VS-DN1": ("DN1_Mean", "DN2_Mean"),
             "Group2_DN-VS-NC": ("NC_Mean", "DN_Mean")}
ENR_MAP = {"Group1_DN1-VS-NC": "Group1中DN1比NC",
           "Group1_DN2-VS-DN1": "Group1中DN2比DN1",
           "Group1_DN2-VS-NC": "Group1中DN2比NC",
           "Group2_DN-VS-NC": "Group2中DN比NC"}

def write(df, name, sources):
    with open(SD / name, "w") as f:
        f.write(f"# source: {sources} | extractor: extract_cosmx.py\n")
        df.to_csv(f, index=False)

def dot2space(s):
    """'Cortical.fibroblasts' -> 'Cortical fibroblasts'; 'A..B' -> 'A -> B'; 'A..B.L..R' -> 'A -> B L -> R'"""
    return s.replace("..", " -> ").replace(".", " ")

def norm_de_list(f, entity_col, comp):
    """Normalize a groups/ DE-list CSV to fixed schema."""
    raw = pd.read_csv(f)
    raw = raw.rename(columns={raw.columns[0]: "entity"})
    ctrl_col, case_col = MEAN_COLS[comp]
    out = pd.DataFrame({
        "comparison": comp, entity_col: raw["entity"],
        "W": raw["W"], "pvalue": raw["Pvalue"], "significance": raw["Significance"],
        "NC_mean": raw[ctrl_col], "case_mean": raw[case_col],
        "mean_diff": raw["Mean_difference"], "change": raw["Change"]})
    out[entity_col] = out[entity_col].map(dot2space)
    return out

def main():
    # 1) samples + groups
    g = pd.read_csv(SRC / "groups/groups.csv", index_col=0).reset_index()
    g = g.rename(columns={g.columns[0]: "sample"})
    md = pd.read_csv(SRC / "gene_NMF/cNMF.metadata.csv")
    md["zone"] = md["zone_leiden"].apply(lambda z: f"Zone {int(z)}" if pd.notna(z) else "Zone -1")
    n_cells = md.groupby("sample").size().rename("n_cells")
    fovs = pd.read_csv(SRC / "YKKY0385_20241121_fov_sample.csv").groupby("sample")["fov"].nunique().rename("n_fovs")
    samples = g[["sample", "Group1", "Group2"]].rename(columns={"Group1": "group1", "Group2": "group2"})
    samples = samples.merge(n_cells, on="sample", how="left").merge(fovs, on="sample", how="left")
    samples["n_cells"] = samples["n_cells"].astype("Int64"); samples["n_fovs"] = samples["n_fovs"].astype("Int64")
    write(samples, "samples.csv", "groups/groups.csv; YKKY0385_20241121_fov_sample.csv; cNMF.metadata.csv")

    # 2) cells.csv (metadata: 18 types, no Mix) + enrich_zscores
    cells = md[["sample", "x_slide_mm", "y_slide_mm", "main_celltype", "new_celltype", "zone", "fov"]].copy()
    cells["major_group"] = cells["main_celltype"].map(GROUPS_MAP)
    cells = cells.merge(samples[["sample", "group1"]], on="sample")
    cells = cells.rename(columns={"x_slide_mm": "x_mm", "y_slide_mm": "y_mm"})
    write(cells, "cells.csv", "cNMF.metadata.csv")
    enr = pd.read_csv(SRC / "cell_colocate/celltype_enrich/allsample.enrich.csv", index_col=0).reset_index()
    enr = enr.rename(columns={enr.columns[0]: "sample"})
    enr.columns = [dot2space(c) if "->" not in c else c for c in enr.columns]
    enr = enr.merge(samples[["sample", "group1"]], on="sample", how="left")
    write(enr, "enrich_zscores.csv", "cell_colocate/celltype_enrich/allsample.enrich.csv")

    # 3) h5ad: umap + dotplot + gene-module expression (chunked, backed mode)
    import anndata as ad
    a = ad.read_h5ad(SRC / "cell_anno/filter.cell_anno.h5ad", backed="r")
    obs = a.obs[["sample", "main_celltype", "new_celltype"]]
    um = pd.DataFrame(a.obsm["banksy_umap"], columns=["UMAP1", "UMAP2"])
    um["sample"] = obs["sample"].values
    um["main_celltype"] = obs["main_celltype"].map(GROUPS_MAP)
    um["new_celltype"] = obs["new_celltype"].values
    um["major_group"] = um["main_celltype"]
    write(um, "umap.csv", "filter.cell_anno.h5ad (banksy_umap; 53,883 cells incl Mix)")

    genes = ["PODXL", "WT1", "PECAM1", "CDH5", "SLC12A3", "AQP3", "COL1A1", "PDGFRA",
             "CD3D", "CD4", "CD14", "CD68", "CD79A", "IGKC"]
    vn = list(a.var_names)
    gidx = np.array([vn.index(gg) for gg in genes])
    ct = np.asarray(obs["new_celltype"].values)
    sums = np.zeros((19, len(genes))); nz = np.zeros_like(sums); cnt = np.zeros(19)
    labs = sorted(obs["new_celltype"].unique())
    for i in range(0, a.n_obs, 10000):
        X = a[i:i + 10000, gidx].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        c = ct[i:i + 10000]
        for k, lab in enumerate(labs):
            m = c == lab
            if m.any():
                sums[k] += X[m].sum(0); nz[k] += (X[m] > 0).sum(0); cnt[k] += m.sum()
    lab2group = (obs.groupby("new_celltype")["main_celltype"]
                 .agg(lambda s: s.dropna().iloc[0] if len(s.dropna()) else None)).map(GROUPS_MAP)
    rows = []
    for k, lab in enumerate(labs):
        mean = sums[k] / max(cnt[k], 1); pct = 100 * nz[k] / max(cnt[k], 1)
        for gi, gg in enumerate(genes):
            rows.append(dict(celltype=lab, main_group=lab2group.get(lab, ""),
                             gene=gg, pct_expressed=round(pct[gi], 2), avg_expr_scaled=round(float(mean[gi]), 4)))
    write(pd.DataFrame(rows), "dotplot_values.csv", "filter.cell_anno.h5ad (chunked grouped mean)")

    # 4) proportions
    prop = pd.read_csv(SRC / "cell_anno/sample_celltype_proportion.csv")
    prop = prop.merge(samples[["sample", "group1"]], on="sample", how="left")
    cols = [c for c in common.CELLTYPES if c in prop.columns]
    assert len(cols) == 18, f"proportion columns missing: {set(common.CELLTYPES) - set(cols)}"
    write(prop[["sample", "group1"] + cols], "celltype_proportions.csv", "cell_anno/sample_celltype_proportion.csv")
    mg = prop[["sample", "group1"]].copy()
    for grp in common.GROUP_COLORS:
        member = [c for c in cols if common.MAJOR_GROUP[c] == grp]
        mg[grp] = prop[member].sum(1)
    write(mg, "major_group_proportions.csv", "aggregated from sample_celltype_proportion.csv")

    # 5) zones — zone_markergene.txt is tab-separated: [Label, Higher_exp_markers(comma-separated)]
    zraw = pd.read_csv(SRC / "cell_zone/zone_markergene.txt", sep="\t")
    zrows = []
    for _, r in zraw.iterrows():
        zone = str(r["Label"]).replace("zone_", "Zone ")
        genes_ = [x.strip() for x in str(r["Higher_exp_markers"]).split(",") if x.strip()][:10]
        for gi, gg in enumerate(genes_, 1):
            zrows.append(dict(zone=zone, gene=gg, rank=gi))
    zmark = pd.DataFrame(zrows)
    write(zmark, "zone_markers.csv", "cell_zone/zone_markergene.txt (top 10 per zone)")
    zp = pd.read_csv(SRC / "cell_zone/sample_zone_proportion.csv")
    zp = zp.merge(samples[["sample", "group1"]], on="sample", how="left")
    zcols = {c: f"Zone {int(c)}" for c in zp.columns if str(c).isdigit()}
    zp = zp.rename(columns=zcols)
    write(zp, "zone_proportions.csv", "cell_zone/sample_zone_proportion.csv")
    zc = cells.groupby(["zone", "new_celltype"]).size().rename("n_cells").reset_index()
    write(zc, "zone_composition.csv", "cNMF.metadata.csv")

    # 6) group DE lists (common schema)
    ct_de = pd.concat([norm_de_list(SRC / f"groups/cell_type/proportion/{c}/CellType_prop_DE_list_in_{c}.csv", "celltype", c) for c in COMP_DIRS])
    write(ct_de, "ct_prop_de.csv", "groups/cell_type/proportion/*/CellType_prop_DE_list_in_*.csv")
    enr_de = pd.concat([norm_de_list(SRC / f"groups/cell_colocate/celltype_enrich/{c}/CellEnrich_DE_list_in_{c}.csv", "cell_pair", c) for c in COMP_DIRS])
    write(enr_de, "enrich_de.csv", "groups/cell_colocate/celltype_enrich/*/CellEnrich_DE_list_in_*.csv")
    lr_num = pd.concat([norm_de_list(SRC / f"groups/cell_communication/cellchat/{c}/cell_pair/cell_lr_num_DE_list_in_{c}.csv", "cell_pair", c) for c in COMP_DIRS])
    write(lr_num, "lr_num_de.csv", "groups/cell_communication/cellchat/*/cell_pair/cell_lr_num_DE_list_in_*.csv")
    lr_val = pd.concat([norm_de_list(SRC / f"groups/cell_communication/cellchat/{c}/lr_pair/lr_value_DE_list_in_{c}.csv", "lr_pair", c) for c in COMP_DIRS])
    write(lr_val, "lr_value_de.csv", "groups/cell_communication/cellchat/*/lr_pair/lr_value_DE_list_in_*.csv")
    # 6b) per-sample LR-pair communication probabilities (Fig 6e single-pair panel)
    lr_samp_frames = []
    for c in COMP_DIRS:
        raw = pd.read_csv(SRC / f"groups/cell_communication/cellchat/{c}/lr_pair/lr_value_DE_list_in_{c}.csv")
        raw = raw.rename(columns={raw.columns[0]: "lr_pair"})
        samp_cols = [s for s in raw.columns if s in set(samples["sample"])]
        long = raw.melt(id_vars=["lr_pair"], value_vars=samp_cols,
                        var_name="sample", value_name="value")
        long["comparison"] = c
        long["lr_pair"] = long["lr_pair"].map(dot2space)
        long = long.merge(samples[["sample", "group1"]], on="sample", how="left")
        lr_samp_frames.append(long)
    write(pd.concat(lr_samp_frames, ignore_index=True), "lr_pair_samples.csv",
          "groups/cell_communication/cellchat/*/lr_pair/lr_value_DE_list_in_*.csv (per-sample columns)")
    cent = pd.concat([
        norm_de_list(SRC / f"groups/spatial_distr/centrality/{c}/Cell{m}_DE_list_in_{c}.csv", "celltype", c).assign(metric=m)
        for c in COMP_DIRS for m in ["Degree", "Closeness", "Clustering"]])
    write(cent, "centrality_de.csv", "groups/spatial_distr/centrality/*/Cell*_DE_list_in_*.csv")
    mp_de = pd.concat([norm_de_list(SRC / f"groups/gene_NMF/{c}/cNMF_DE_list_in_{c}.csv", "mp", c) for c in COMP_DIRS])
    write(mp_de, "mp_de.csv", "groups/gene_NMF/*/cNMF_DE_list_in_*.csv")
    zone_de = pd.concat([norm_de_list(SRC / f"groups/cell_zone/proportion/{c}/CellZone_prop_DE_list_in_{c}.csv", "zone", c) for c in COMP_DIRS])
    zone_de["zone"] = zone_de["zone"].str.replace("X", "Zone ")
    write(zone_de, "zone_prop_de.csv", "groups/cell_zone/proportion/*/CellZone_prop_DE_list_in_*.csv")

    # 7) per-celltype edgeR DEG (4 comparisons) -> ct_deg.csv
    frames = []
    for d in sorted(gb.glob(str(SRC / "groups/cell_type/DEG_result/Matrix/*"))):
        if not (Path(d)).is_dir():
            continue
        ct = Path(d).name.replace("_", " ")
        for comp, tag in XLSX_MAP.items():
            f = Path(d) / f"DE_list_{tag}_edgeR.xlsx"
            if f.exists():
                dd = pd.read_excel(f)
                dd = dd.rename(columns={"GeneSymbol": "gene", "logFC": "logFC", "PValue": "PValue", "FDR": "fdr"})
                dd["celltype"] = ct; dd["comparison"] = comp
                frames.append(dd[["celltype", "comparison", "gene", "logFC", "fdr"]])
    ct_deg = pd.concat(frames, ignore_index=True)
    write(ct_deg, "ct_deg.csv", "groups/cell_type/DEG_result/Matrix/*/DE_list_*_edgeR.xlsx")

    # 8) zone DEG + pathways -> zone_deg.csv, zone_pathways.csv
    zframes = []
    for d in sorted(gb.glob(str(SRC / "groups/cell_zone/DEG_result/Matrix/*"))):
        if not Path(d).is_dir():
            continue
        zone = f"Zone {Path(d).name}"
        for comp, tag in XLSX_MAP.items():
            f = Path(d) / f"DE_list_{tag}_edgeR.xlsx"
            if f.exists():
                dd = pd.read_excel(f)
                dd = dd.rename(columns={"GeneSymbol": "gene", "logFC": "logFC", "PValue": "PValue", "FDR": "fdr"})
                dd["zone"] = zone; dd["comparison"] = comp
                zframes.append(dd[["zone", "comparison", "gene", "logFC", "fdr"]])
    write(pd.concat(zframes, ignore_index=True), "zone_deg.csv", "groups/cell_zone/DEG_result/Matrix/*/DE_list_*_edgeR.xlsx")
    pframes = []
    for z in range(11):
        for comp, tag in ENR_MAP.items():
            f = SRC / f"groups/cell_zone/Enrichment_results/{z}/GSEA_KEGG/{tag}_GSEA_KEGG_enrichment.xlsx"
            if f.exists():
                dd = pd.read_excel(f)
                dd["zone"] = f"Zone {z}"; dd["comparison"] = comp
                pframes.append(dd[["zone", "comparison", "ID", "Description", "setSize", "NES", "pvalue", "p.adjust"]])
    write(pd.concat(pframes, ignore_index=True), "zone_pathways.csv", "groups/cell_zone/Enrichment_results/*/GSEA_KEGG/*.xlsx")

    # 9) communication strength (long) from per-sample CC_matrix.csv
    cc = []
    for s in samples["sample"]:
        m = pd.read_csv(SRC / f"cell_communication/cellchat/{s}/{s}.CC_matrix.csv", index_col=0)
        m.index = [dot2space(str(i)) for i in m.index]
        m.columns = [dot2space(str(c)) for c in m.columns]
        mm = m.stack().rename("strength").reset_index()
        mm.columns = ["source", "target", "strength"]
        mm["sample"] = s
        cc.append(mm)
    cc = pd.concat(cc, ignore_index=True).merge(samples[["sample", "group1"]], on="sample")
    write(cc, "comm_strength.csv", "cell_communication/cellchat/*/CC_matrix.csv (long)")

    # 10) NMF
    mp = pd.read_csv(SRC / "gene_NMF/mean_MPs_sample.csv").drop(columns=["MP4"], errors="ignore")
    mp = mp.merge(samples[["sample", "group1"]], on="sample", how="left")
    write(mp, "mp_scores.csv", "gene_NMF/mean_MPs_sample.csv (MP4 dropped)")
    write(pd.read_csv(SRC / "gene_NMF/MP_gene.csv"), "mp_genes.csv", "gene_NMF/MP_gene.csv")
    write(pd.read_csv(SRC / "gene_NMF/MP_C2_pathway.csv"), "mp_pathways.csv", "gene_NMF/MP_C2_pathway.csv")

    # 11) gene modules for Fig4B: union of top-30 |logFC| genes per celltype (Group2 DN-vs-NC),
    #     per-sample mean expression (chunked over h5ad), z-scored per gene, ward k=3
    top_genes = (ct_deg[ct_deg.comparison == "Group2_DN-VS-NC"]
                 .assign(a=lambda d: d.logFC.abs())
                 .sort_values(["fdr", "a"], ascending=[True, False])
                 .groupby("celltype").head(30)["gene"].unique())
    top_genes = [gg for gg in top_genes if gg in vn]
    tg_idx = np.array([vn.index(gg) for gg in top_genes])
    sarr = np.asarray(obs["sample"].values)
    sample_list = samples["sample"].tolist()
    acc = {s: np.zeros(len(top_genes)) for s in sample_list}
    for i in range(0, a.n_obs, 10000):
        X = a[i:i + 10000, tg_idx].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        ss = sarr[i:i + 10000]
        for j, s in enumerate(sample_list):
            m = ss == s
            if m.any():
                acc[s] += X[m].sum(0)
    del a
    gsm = pd.DataFrame(acc, index=top_genes).T
    gsm.index.name = "sample"
    zs = (gsm - gsm.mean(0)) / gsm.std(0)
    zlong = zs.reset_index().melt(id_vars="sample", var_name="gene", value_name="zscore")
    zlong = zlong.merge(samples[["sample", "group1"]], on="sample")
    write(zlong, "gene_module_zscores.csv", "filter.cell_anno.h5ad + ct_deg.csv (top-30/celltype union)")
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist
    Z = linkage(zs.T, method="ward")
    cl = fcluster(Z, t=3, criterion="maxclust")
    assign = pd.DataFrame({"gene": top_genes, "module_id": cl})
    # name modules by content: draft expects up-module to contain CD74/COL1A2/DCN/LYZ/C1S
    content = {k: set(g) for k, g in assign.groupby("module_id")["gene"]}
    def name_module(gset):
        if {"CD74", "COL1A2", "DCN", "LYZ", "C1S"} & gset:
            return "up"
        if gset & {"PCK1", "GATM", "APOE", "GPX3", "ALDOB", "SLC13A3", "ASS1", "FABP1"}:
            return "down"
        return "mixed"
    assign["module"] = assign["module_id"].map(lambda k: name_module(content[k]))
    write(assign[["gene", "module"]], "gene_module_assignment.csv", "derived (ward k=3 on z-scored sample means)")

    # 12) QC
    qc = pd.read_csv(SRC / "QC/YKKY0385_20241121.qc_overview.xls", sep="\t", header=None, names=["metric", "value"])
    write(qc, "qc_overview.csv", "QC/YKKY0385_20241121.qc_overview.xls")
    qf = pd.read_csv(SRC / "QC/YKKY0385_20241121.qc_fovs.xls", sep="\t")
    write(qf, "qc_fovs.csv", "QC/YKKY0385_20241121.qc_fovs.xls")
    fq = pd.read_excel(SRC / "QC/filter_qc_df.xlsx").iloc[0]
    fq_long = pd.DataFrame({"metric": [str(i) for i in fq.index[1:]],
                            "value": [float(v) for v in fq.values[1:]]})
    write(fq_long, "filter_qc.csv", "QC/filter_qc_df.xlsx")
    print("done:", len(list(SD.glob("*.csv"))), "files")

if __name__ == "__main__":
    main()
