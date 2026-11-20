"""Fig 4: group-level remodelling (183 mm). Stars-only heatmaps, dedicated legends, clean boxplot."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import common
from common import SD, FIGDIR

np.random.seed(42)
common.set_nc_style()
ct_de = pd.read_csv(SD / "ct_prop_de.csv", comment="#")
enr_de = pd.read_csv(SD / "enrich_de.csv", comment="#")
zs = pd.read_csv(SD / "gene_module_zscores.csv", comment="#")
ass = pd.read_csv(SD / "gene_module_assignment.csv", comment="#")
enr_z = pd.read_csv(SD / "enrich_zscores.csv", comment="#")

COMPS = ["Group1_DN1-VS-NC", "Group1_DN2-VS-NC", "Group2_DN-VS-NC"]
COMP_LAB = {"Group1_DN1-VS-NC": "DN1 vs NC", "Group1_DN2-VS-NC": "DN2 vs NC", "Group2_DN-VS-NC": "DN vs NC"}
GRP_COL = {"NC": "#999999", "DN1": "#0072B2", "DN2": "#D55E00"}
MOD_COL = {"up": "#B2182B", "down": "#2166AC"}

fig = plt.figure(figsize=(common.W_DOUBLE, 7.4))
fig.subplots_adjust(bottom=0.16)   # room for rotated 90deg pair labels in panel c
gs = fig.add_gridspec(2, 2, hspace=0.72, wspace=0.45, height_ratios=[1, 1.1])

def star(v, s, vmax):
    return dict(text="*" if s == "Sig" else "", color="white" if abs(v) > vmax * 0.55 else "black")

# --- a: celltype proportion DE heatmap
axA = fig.add_subplot(gs[0, 0]); common.panel_label(axA, "a")
dd = ct_de[ct_de.comparison.isin(COMPS)].copy()
dd["cmp"] = dd.comparison.map(COMP_LAB)
mat = dd.pivot_table(index="cmp", columns="celltype", values="mean_diff", aggfunc="first").reindex(columns=common.CELLTYPES, index=list(COMP_LAB.values()))
smat = dd.pivot_table(index="cmp", columns="celltype", values="significance", aggfunc="first").reindex(columns=common.CELLTYPES, index=list(COMP_LAB.values()))
vmax = np.nanmax(np.abs(mat.values)) or 1
im = axA.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
axA.set_xticks(range(18)); axA.set_xticklabels([common.CT_SHORT[c] for c in common.CELLTYPES], rotation=90, fontsize=6.5)
axA.set_yticks(range(3)); axA.set_yticklabels(mat.index, fontsize=7.5)
for i in range(3):
    for j in range(18):
        v = mat.values[i, j]
        if pd.notna(v):
            axA.text(j, i, star(v, smat.values[i, j], vmax)["text"], ha="center", va="center",
                     fontsize=8, color=star(v, smat.values[i, j], vmax)["color"])
cb = plt.colorbar(im, ax=axA, shrink=0.75, pad=0.02); cb.set_label("Δ proportion (%)", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)

# --- b: gene module heatmap
axB = fig.add_subplot(gs[0, 1]); common.panel_label(axB, "b")
piv = zs.pivot_table(index="gene", columns="sample", values="zscore")
piv["mod"] = piv.index.map(lambda g: ass.set_index("gene").module.get(g, "down"))
mod_rank = {"up": 0, "down": 1}
piv = piv.sort_values("mod", key=lambda s: s.map(mod_rank))
mod_row = piv["mod"].map(MOD_COL)
samples = [c for c in piv.columns if c != "mod"]
grp_cols = [GRP_COL[zs[zs["sample"] == s].group1.iloc[0]] for s in samples]
pivm = piv[samples]
vmax2 = np.nanmax(np.abs(pivm.values)) or 1
im2 = axB.imshow(pivm.values, cmap="RdBu_r", vmin=-vmax2, vmax=vmax2, aspect="auto", interpolation="nearest")
axB.set_xticks(range(len(samples))); axB.set_xticklabels(samples, fontsize=6.8, rotation=90)
for i, c in enumerate(grp_cols):
    axB.axvspan(i - 0.5, i + 0.5, ymin=0.0, ymax=0.022, color=c, lw=0)
axB.set_yticks([]); axB.set_ylabel("Genes (top-30 DEGs per cell type)", fontsize=7.5)
cb2 = plt.colorbar(im2, ax=axB, shrink=0.75, pad=0.02); cb2.set_label("z-scored mean expr.", fontsize=7.5); cb2.ax.tick_params(labelsize=6.5)
mh = [mpatches.Patch(fc=c, label=m) for m, c in MOD_COL.items()]
common.hline_legend(axB, mh, None, ncol=2)

# --- c: colocalization enrichment DE heatmap (top 25 pairs)
axC = fig.add_subplot(gs[1, 0]); common.panel_label(axC, "c")
ee = enr_de[enr_de.comparison.isin(COMPS)].copy()
ee["cmp"] = ee.comparison.map(COMP_LAB)
top_pairs = (ee.groupby("cell_pair")["mean_diff"].agg(lambda s: s.abs().max())
             .sort_values(ascending=False).head(25).index.tolist())
mat3 = ee.pivot_table(index="cmp", columns="cell_pair", values="mean_diff", aggfunc="first").reindex(index=list(COMP_LAB.values()), columns=top_pairs)
smat3 = ee.pivot_table(index="cmp", columns="cell_pair", values="significance", aggfunc="first").reindex(index=list(COMP_LAB.values()), columns=top_pairs)
vmax3 = np.nanmax(np.abs(mat3.values)) or 1
im3 = axC.imshow(mat3.values, cmap="RdBu_r", vmin=-vmax3, vmax=vmax3, aspect="auto")
def short_pair(p):
    """'SRC->TGT' with nephrology abbreviations, e.g. 'M2 macro→Monocyte'."""
    a, b = [x.strip() for x in p.split("->")]
    return f"{common.CT_ABBR[a]}→{common.CT_ABBR[b]}"

axC.set_xticks(range(25)); axC.set_xticklabels([short_pair(p) for p in top_pairs], rotation=90, fontsize=6)
axC.set_yticks(range(3)); axC.set_yticklabels(mat3.index, fontsize=7.5)
for i in range(3):
    for j in range(25):
        v = mat3.values[i, j]
        if pd.notna(v):
            axC.text(j, i, star(v, smat3.values[i, j], vmax3)["text"], ha="center", va="center",
                     fontsize=7, color=star(v, smat3.values[i, j], vmax3)["color"])
cb3 = plt.colorbar(im3, ax=axC, shrink=0.75, pad=0.02); cb3.set_label("Δ enrichment z", fontsize=7.5); cb3.ax.tick_params(labelsize=6.5)

# --- d: representative pair boxplot
axD = fig.add_subplot(gs[1, 1]); common.panel_label(axD, "d")
PAIR = "M2 Macrophage->Monocyte cell"
assert PAIR in enr_z.columns, f"pair column missing: {PAIR}"
bp = enr_z[["sample", "group1", PAIR]].dropna()
order = ["NC", "DN1", "DN2"]
data = [bp[bp.group1 == g][PAIR].values for g in order]
bp_art = axD.boxplot(data, positions=[0, 1, 2], widths=0.45, patch_artist=True,
                     medianprops=dict(color="black", lw=1), flierprops=dict(markersize=2))
for patch, g in zip(bp_art["boxes"], order):
    patch.set_facecolor(GRP_COL[g]); patch.set_alpha(0.85)
for xi, g in enumerate(order):
    yy = bp[bp.group1 == g][PAIR].values
    axD.scatter(np.random.normal(xi, 0.05, len(yy)), yy, s=5, c="#444444", alpha=0.6, zorder=3, edgecolors="none")
axD.set_xticks([0, 1, 2]); axD.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(order, data)], fontsize=7.5)
axD.set_ylabel("Colocalization enrichment (z)")
axD.set_title("M2 macrophage → monocyte", fontsize=8)

common.save_fig(fig, "Fig4_group_remodelling")
print("saved Fig4_group_remodelling")
