"""Fig 1: single-cell and spatial atlas of the DKD kidney (183 mm).
Layout: row1 a(schematic)|b(UMAP 8 groups)|c(UMAP 18 types)|legend-column;
row2 d(dotplot)|e(spatial rep + 16-sample grid); row3 f(proportions)."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import common
from common import SD, FIGDIR

common.set_nc_style()
df_umap = pd.read_csv(SD / "umap.csv", comment="#")
df_cells = pd.read_csv(SD / "cells.csv", comment="#")
df_dot = pd.read_csv(SD / "dotplot_values.csv", comment="#")
df_mg = pd.read_csv(SD / "major_group_proportions.csv", comment="#").set_index("sample")
df_samples = pd.read_csv(SD / "samples.csv", comment="#")

MIX_COLOR = "#BBBBBB"
DOT = dict(linewidths=0)

fig = plt.figure(figsize=(common.W_DOUBLE, 9.7))   # 246 mm < 247 mm full-page limit
gs = fig.add_gridspec(3, 4, height_ratios=[1, 1.05, 0.8], hspace=0.5, wspace=0.42,
                      width_ratios=[0.72, 0.94, 0.94, 0.90])

# --- a: study design (vertical flow — fits the narrow column without overlaps)
axA = fig.add_subplot(gs[0, 0]); common.panel_label(axA, "a")
axA.axis("off")
steps = [("Tissue microarray", "16 samples\n5 DKD + 11 control"),
         ("CosMx SMI 1,000-plex", "47 fields of view"),
         ("QC & clustering", "61,998 → 53,883 cells"),
         ("Atlas", "18 cell types\n8 groups · 11 zones")]
axA.set_xlim(0, 1); axA.set_ylim(0, 1)
for i, (head, sub) in enumerate(steps):
    ytop = 0.97 - i * 0.25
    axA.annotate("", xy=(0.5, ytop - 0.085), xytext=(0.5, ytop - 0.028),
                 arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.0),
                 annotation_clip=False)
    axA.text(0.5, ytop, head, ha="center", va="top", fontsize=7.5, fontweight="bold")
    axA.text(0.5, ytop - 0.026, sub, ha="center", va="top", fontsize=7, color="#444444")

# --- b: UMAP by 8 major groups (legend below, dedicated space)
axB = fig.add_subplot(gs[0, 1]); common.panel_label(axB, "b")
handles = []
for grp, col in common.GROUP_COLORS.items():
    m = df_umap.major_group == grp
    axB.scatter(df_umap.loc[m, "UMAP1"], df_umap.loc[m, "UMAP2"], s=1.4, c=col, rasterized=True, **DOT)
    handles.append(mpatches.Patch(fc=col, label=grp))
axB.set(xlabel="UMAP1", ylabel="UMAP2", xticks=[], yticks=[])
common.hline_legend(axB, handles, None, ncol=2)

# --- c: UMAP by 18 cell types
axC = fig.add_subplot(gs[0, 2]); common.panel_label(axC, "c")
for ct, col in common.CELLTYPE_COLORS.items():
    m = df_umap.new_celltype == ct
    axC.scatter(df_umap.loc[m, "UMAP1"], df_umap.loc[m, "UMAP2"], s=1.4, c=col, rasterized=True, **DOT)
m_mix = df_umap.new_celltype == "Mix"
axC.scatter(df_umap.loc[m_mix, "UMAP1"], df_umap.loc[m_mix, "UMAP2"], s=1.4, c=MIX_COLOR, rasterized=True, **DOT)
axC.set(xlabel="UMAP1", ylabel="UMAP2", xticks=[], yticks=[])

# --- legend column (dedicated space, single column of 19 rows, abbreviated names)
axL = fig.add_subplot(gs[0, 3]); axL.axis("off")
axL.set_xlim(0, 1); axL.set_ylim(0, 1)
items = list(common.CELLTYPE_COLORS.items()) + [("Mix", MIX_COLOR)]
pitch = 1.0 / len(items)
for i, (ct, col) in enumerate(items):
    yc = 0.975 - (i + 0.5) * pitch
    axL.add_patch(mpatches.FancyBboxPatch((0.02, yc - 0.022), 0.16, 0.044,
                 boxstyle="round,pad=0.008", fc=col, ec="none"))
    axL.text(0.22, yc, common.CT_SHORT[ct] if ct != "Mix" else "Mix", fontsize=6.5, va="center")
axL.set_title("Cell types", fontsize=8, loc="left")

# --- d: marker dotplot (14 genes x 18 types)
axD = fig.add_subplot(gs[1, 0:2]); common.panel_label(axD, "d")
dd = df_dot[df_dot.celltype.isin(common.CELLTYPES)]
piv_pct = dd.pivot(index="gene", columns="celltype", values="pct_expressed")
piv_avg = dd.pivot(index="gene", columns="celltype", values="avg_expr_scaled")
order = list(common.CELLTYPE_COLORS)
piv_pct, piv_avg = piv_pct[order], piv_avg[order]
ymax = piv_avg.max().max() or 1
sc = axD.scatter(np.repeat(np.arange(18), piv_pct.shape[0]),
                 np.tile(np.arange(piv_pct.shape[0]), 18),
                 s=piv_pct.T.values.ravel() * 1.6,
                 c=piv_avg.T.values.ravel(), cmap="viridis", vmin=0, vmax=ymax,
                 edgecolors="none")
axD.set_xticks(range(18)); axD.set_xticklabels([common.CT_SHORT[c] for c in order],
                                               rotation=90, fontsize=6.8)
axD.set_yticks(range(piv_pct.shape[0])); axD.set_yticklabels(piv_pct.index, fontsize=7.5)
cb = plt.colorbar(sc, ax=axD, shrink=0.75, pad=0.02); cb.set_label("Mean expr.", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)

# --- e: spatial maps — representative large + 16-sample grid
gsE = gs[1, 2:].subgridspec(1, 2, wspace=0.3)
axE = fig.add_subplot(gsE[0]); common.panel_label(axE, "e")
rep = df_cells.groupby("sample").size().sort_values().index[8]
sub = df_cells[df_cells["sample"] == rep]
for grp, col in common.GROUP_COLORS.items():
    m = sub.major_group == grp
    axE.scatter(sub.loc[m, "x_mm"], sub.loc[m, "y_mm"], s=2, c=col, rasterized=True, **DOT)
axE.set_title(f"{rep} (representative)", fontsize=8)
axE.set(xlabel="x (mm)", ylabel="y (mm)"); axE.set_aspect("equal"); axE.tick_params(labelsize=7)
axE.xaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))   # no out-of-view ticks
axE.yaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))

gsM = gsE[1].subgridspec(4, 4, hspace=0.32, wspace=0.08)
for i, s in enumerate(df_samples["sample"].tolist()):
    axm = fig.add_subplot(gsM[i // 4, i % 4])
    ssub = df_cells[df_cells["sample"] == s]
    for grp, col in common.GROUP_COLORS.items():
        m = ssub.major_group == grp
        axm.scatter(ssub.loc[m, "x_mm"], ssub.loc[m, "y_mm"], s=0.45, c=col, rasterized=True, **DOT)
    axm.set_title(s, fontsize=6.5, pad=1.5); axm.set_aspect("equal")
    axm.set_xticks([]); axm.set_yticks([])

# --- f: stacked proportions per sample
axF = fig.add_subplot(gs[2, :]); common.panel_label(axF, "f")
order_s = df_mg.index.tolist()
bottom = np.zeros(16)
handles = []
for grp, col in common.GROUP_COLORS.items():
    axF.bar(range(16), df_mg[grp].values, bottom=bottom, color=col, width=0.72, label=grp)
    bottom += df_mg[grp].values
    handles.append(mpatches.Patch(fc=col, label=grp))
axF.set_xticks(range(16)); axF.set_xticklabels(order_s, fontsize=7.5)
axF.set_ylabel("% cells"); axF.set_ylim(0, 100)
common.hline_legend(axF, handles, None, ncol=8)

common.save_fig(fig, "Fig1_atlas")
print("saved Fig1_atlas")
