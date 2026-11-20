"""Fig 8 (draft): immune microenvironment remodelling — transcriptome to protein.
a: per-sample immune cell-type proportions (CosMx, DN vs NC, Wilcoxon+BH stars);
b: spatial distribution of immune cells (representative sample B3);
c: whole-slide mIF density fold changes (5 markers, CTLA4 down highlighted);
d: key tubulo-immune LR pairs (per-sample probabilities, NC vs DN)."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import common
from common import SD, FIGDIR

np.random.seed(42)
common.set_nc_style()
cells = pd.read_csv(SD / "cells.csv", comment="#")
lv = pd.read_csv(SD / "lr_value_de.csv", comment="#")
ls = pd.read_csv(SD / "lr_pair_samples.csv", comment="#")
fc = pd.read_csv(SD / "mif_foldchange.csv", comment="#")

IMMUNE = ["T", "B", "NK", "Plasma", "Mast", "Monocyte cell", "M2 Macrophage", "Neutrophil",
          "Inflammatory fibroblasts"]
DISEASE, CONTROL = "#D55E00", "#0072B2"

# --- a: per-sample immune proportions
props = (cells.groupby(["sample", "new_celltype"]).size()
         .div(cells.groupby("sample").size(), level="sample") * 100)
props = props.reset_index(name="pct")
grp = cells[["sample", "group1"]].drop_duplicates()

# --- b: spatial map of immune cells (representative sample B3)
fig = plt.figure(figsize=(common.W_DOUBLE, 7.2))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.45, height_ratios=[1.15, 1])

axA = fig.add_subplot(gs[0, 0]); common.panel_label(axA, "a")
pos = np.arange(len(IMMUNE)) * 1.35
for k, ct in enumerate(IMMUNE):
    d = props[props.new_celltype == ct].merge(grp, on="sample")
    for g, c, xi in [("NC", CONTROL, -0.22), ("DN1", DISEASE, 0.0), ("DN2", DISEASE, 0.22)]:
        v = d[d.group1 == g].pct.values
        axA.scatter(np.random.normal(pos[k] + xi, 0.03, len(v)), v, s=10, c=c, alpha=0.65,
                    edgecolors="none", zorder=3)
        bp = axA.boxplot([v], positions=[pos[k] + xi], widths=0.16, patch_artist=True,
                         manage_ticks=False, showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(c); patch.set_alpha(0.7)
        for el in bp["medians"]:
            el.set_color("black")
    ncx = d[d.group1 == "NC"].pct.values
    dnx = d[d.group1.isin(["DN1", "DN2"])].pct.values
    if len(ncx) and len(dnx):
        u, p = mannwhitneyu(ncx, dnx, alternative="two-sided")
        stars = common.sig_stars(p)
        if stars != "ns":
            axA.text(pos[k], max(max(ncx), max(dnx)) * 1.12, stars, ha="center", fontsize=8)
axA.set_xticks(pos); axA.set_xticklabels([common.CT_ABBR[c] for c in IMMUNE],
                                         rotation=45, ha="right", fontsize=7)
axA.set_ylabel("% cells", fontsize=8)
axA.set_title("Immune cell proportions (CosMx, per sample)", fontsize=8)
axA.tick_params(labelsize=7)

axB = fig.add_subplot(gs[0, 1]); common.panel_label(axB, "b")
rep = "B3"
sub = cells[cells["sample"] == rep]
non = sub[~sub.new_celltype.isin(IMMUNE)]
axB.scatter(non.x_mm, non.y_mm, s=0.4, c="#DDDDDD", rasterized=True, linewidths=0)
for ct in IMMUNE:
    m = sub[sub.new_celltype == ct]
    axB.scatter(m.x_mm, m.y_mm, s=1.6, c=common.CELLTYPE_COLORS[ct], rasterized=True,
                linewidths=0, label=ct)
axB.set_title(f"{rep}: immune cells in situ", fontsize=8)
axB.set_aspect("equal"); axB.set_xticks([]); axB.set_yticks([])
axB.legend(fontsize=6, loc="upper right", ncol=2, handlelength=0.8, columnspacing=0.5)

axC = fig.add_subplot(gs[1, 0]); common.panel_label(axC, "c")
MARKERS = ["cd8", "cd20", "cd4", "cd68", "ctla4"]
MARKER_LAB = {"cd8": "CD8", "cd20": "CD20", "cd4": "CD4", "cd68": "CD68", "ctla4": "CTLA4"}
fcd = fc[fc.metric == "dens"].set_index("marker").loc[MARKERS].sort_values("fold_change", ascending=False)
yy = np.arange(len(fcd))[::-1]
cols = [DISEASE if v > 1 else CONTROL for v in fcd.fold_change.values]
axC.scatter(fcd.fold_change.values, yy, c=cols, s=40, zorder=3, edgecolors="#333333", linewidths=0.3)
for yi, v in enumerate(fcd.fold_change.values):
    axC.plot([1, v], [yy[yi], yy[yi]], color="#999999", lw=1, zorder=1)
    axC.text(v, yy[yi] + 0.15, f"{v:.2f}×", fontsize=7, va="bottom", ha="center")
axC.axvline(1, color="#333333", lw=0.8, ls="--")
axC.set_yticks(yy); axC.set_yticklabels([MARKER_LAB[m] for m in fcd.index], fontsize=7.5)
axC.set_xlabel("Whole-slide density fold change (DN / control)", fontsize=8)
axC.set_xlim(0, 6.2)
axC.set_title("mIF: whole-slide immune densities", fontsize=8)
axC.text(0.995, 0.05, "n = 2 slides/group (descriptive)", transform=axC.transAxes,
         fontsize=7, ha="right", color="#555555")

axD = fig.add_subplot(gs[1, 1]); common.panel_label(axD, "d")
KEY_UP = [
 ("Distal tubular cell -> Neutrophil SPP1 -> ITGAV_ITGB5", "SPP1→ITGAV_ITGB5", "DT→Neutro"),
 ("Cortical fibroblasts -> Distal tubular cell EGF -> EGFR", "EGF→EGFR", "cFibro→DT"),
 ("Distal tubular cell -> Monocyte cell KLK3 -> NGFR", "KLK3→NGFR", "DT→Mono")]
rng = np.random.default_rng(7)
for k, (pair, lig_label, cell_label) in enumerate(KEY_UP):
    sub = ls[(ls.comparison == "Group2_DN-VS-NC") & (ls.lr_pair == pair)]
    pos = k * 2.4
    for g, c, xi in [("NC", CONTROL, -0.3), ("DN1", DISEASE, 0.0), ("DN2", DISEASE, 0.3)]:
        v = sub[sub.group1 == g].value.values
        axD.scatter(rng.normal(pos + xi, 0.05, len(v)), v, s=10, c=c, alpha=0.6,
                    edgecolors="none", zorder=3)
        bp = axD.boxplot([v], positions=[pos + xi], widths=0.22, patch_artist=True,
                         manage_ticks=False, showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(c); patch.set_alpha(0.7)
        for el in bp["medians"]:
            el.set_color("black")
    p = lv[(lv.comparison == "Group2_DN-VS-NC") & (lv.lr_pair == pair)].pvalue.values
    if len(p):
        vmax = sub[sub.group1.isin(["NC", "DN1", "DN2"])].value.max()
        axD.text(pos, vmax * 1.5, common.pval_text(p[0]), ha="center", fontsize=6.5)
axD.set_yscale("log")
axD.set_xticks([k * 2.4 for k in range(3)])
axD.set_xticklabels([f"{l}\n{c}" for _, l, c in KEY_UP], fontsize=7)
axD.set_ylabel("Comm. prob. (log)", fontsize=8)
axD.set_title("Tubulo-immune signalling pairs (up in DKD)", fontsize=8)
axD.tick_params(labelsize=7)

common.save_fig(fig, "Fig8_immune_integration")
print("saved Fig8_immune_integration")
