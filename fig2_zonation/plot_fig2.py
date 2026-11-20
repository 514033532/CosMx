"""Fig 2: functional zonation (183 mm). Clean layout: legends in dedicated space,
heatmap cells carry significance stars only (values in source_data)."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import common
from common import SD, FIGDIR

common.set_nc_style()
df_cells = pd.read_csv(SD / "cells.csv", comment="#")
zm = pd.read_csv(SD / "zone_markers.csv", comment="#")
zc = pd.read_csv(SD / "zone_composition.csv", comment="#")
zpd = pd.read_csv(SD / "zone_prop_de.csv", comment="#")
zd = pd.read_csv(SD / "zone_deg.csv", comment="#")

FDR_T, FC_T = 0.05, 0.25
DOT = dict(linewidths=0)

fig = plt.figure(figsize=(common.W_DOUBLE, 7.0))
gs = fig.add_gridspec(2, 3, hspace=0.72, wspace=0.5, width_ratios=[1.15, 1.1, 1.1])

# --- a: representative-sample zone map + legend below
axA = fig.add_subplot(gs[0, 0]); common.panel_label(axA, "a")
rep = df_cells.groupby("sample").size().sort_values().index[8]
sub = df_cells[df_cells["sample"] == rep]
for z, col in common.ZONE_COLORS.items():
    m = sub.zone == z
    axA.scatter(sub.loc[m, "x_mm"], sub.loc[m, "y_mm"], s=1.8, c=col, rasterized=True, **DOT)
axA.set_title(f"{rep} (representative)", fontsize=8)
axA.set(xlabel="x (mm)", ylabel="y (mm)"); axA.set_aspect("equal"); axA.tick_params(labelsize=7)
axA.xaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))   # no out-of-view ticks
axA.yaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))
zh = [mpatches.Patch(fc=c, label=z.replace("Zone ", "")) for z, c in common.ZONE_COLORS.items()]
common.hline_legend(axA, zh, None, ncol=6)

# --- b: zone composition (100% stacked by major group)
axB = fig.add_subplot(gs[0, 1]); common.panel_label(axB, "b")
zc2 = zc.merge(zc.groupby("zone")["n_cells"].sum().rename("tot"), on="zone")
zc2["group"] = zc2.new_celltype.map(common.MAJOR_GROUP)
comp = zc2.groupby(["zone", "group"])["n_cells"].sum().unstack("group").fillna(0)
comp = comp.div(comp.sum(1), axis=0) * 100
order_z = [f"Zone {i}" for i in range(11)]
comp = comp.reindex(order_z)
bottom = np.zeros(11)
gh = []
for grp, col in common.GROUP_COLORS.items():
    if grp in comp.columns:
        axB.bar(range(11), comp[grp].values, bottom=bottom, color=col, width=0.72, label=grp)
        bottom += comp[grp].values
        gh.append(mpatches.Patch(fc=col, label=grp))
axB.set_xticks(range(11)); axB.set_xticklabels([z.replace("Zone ", "") for z in order_z], fontsize=7)
axB.set(xlabel="Zone", ylabel="% cells", ylim=(0, 100))
common.hline_legend(axB, gh, None, ncol=4)

# --- c: zone marker dotplot (top 5 per zone — 54 genes; full right column height)
axC = fig.add_subplot(gs[:, 2]); common.panel_label(axC, "c")
zgenes = zm.pivot_table(index="gene", columns="zone", values="rank", aggfunc="first").reindex(columns=order_z)
keep = sorted({g for z in order_z for g in zm[zm.zone == z].nsmallest(5, "rank").gene})
zgenes = zgenes.loc[zgenes.index.isin(keep)]
ys = np.arange(zgenes.shape[0])
for zi, z in enumerate(order_z):
    for gi, gene in enumerate(zgenes.index):
        r = zgenes.loc[gene, z]
        if pd.notna(r):
            axC.scatter(zi, gi, s=(11 - r) * 5.5 + 7, color=common.ZONE_COLORS[z], edgecolors="none", alpha=0.85)
axC.set_xticks(range(11)); axC.set_xticklabels([z.replace("Zone ", "") for z in order_z], fontsize=7)
axC.set_yticks(ys); axC.set_yticklabels(zgenes.index, fontsize=6.5)
axC.set_xlabel("Zone"); axC.invert_yaxis()

# --- d: zone proportion DE heatmap (stars only)
axD = fig.add_subplot(gs[1, 0]); common.panel_label(axD, "d")
comps = ["Group1_DN1-VS-NC", "Group1_DN2-VS-NC", "Group2_DN-VS-NC"]
dd = zpd[zpd.comparison.isin(comps)].copy()
dd["cmp_short"] = dd.comparison.map({"Group1_DN1-VS-NC": "DN1 vs NC", "Group1_DN2-VS-NC": "DN2 vs NC", "Group2_DN-VS-NC": "DN vs NC"})
mat = dd.pivot_table(index="zone", columns="cmp_short", values="mean_diff", aggfunc="first").reindex(index=order_z, columns=["DN1 vs NC", "DN2 vs NC", "DN vs NC"])
smat = dd.pivot_table(index="zone", columns="cmp_short", values="significance", aggfunc="first").reindex(index=order_z, columns=["DN1 vs NC", "DN2 vs NC", "DN vs NC"])
vmax = np.nanmax(np.abs(mat.values)) or 1
im = axD.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
axD.set_xticks(range(3)); axD.set_xticklabels(mat.columns, fontsize=7.5)
axD.set_yticks(range(11)); axD.set_yticklabels([z.replace("Zone ", "") for z in order_z], fontsize=7.5)
for i in range(3):
    for j in range(11):
        if smat.values[j, i] == "Sig":
            v = mat.values[j, i]
            axD.text(i, j, "*", ha="center", va="center", fontsize=8,
                     color="white" if abs(v) > vmax * 0.55 else "black")
cb = plt.colorbar(im, ax=axD, shrink=0.8, pad=0.02); cb.set_label("Mean diff. (% cells)", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)

# --- e: per-zone DEG counts
axE = fig.add_subplot(gs[1, 1]); common.panel_label(axE, "e")
zde = zd[zd.comparison == "Group2_DN-VS-NC"]
sig = zde[(zde.fdr < FDR_T) & (zde.logFC.abs() > FC_T)]
up = sig[sig.logFC > 0].groupby("zone").size().reindex(order_z).fillna(0)
dn = sig[sig.logFC < 0].groupby("zone").size().reindex(order_z).fillna(0)
xx = np.arange(11)
axE.bar(xx, dn.values, color="#2166AC", width=0.72, label="Down in DN")
axE.bar(xx, up.values, bottom=dn.values, color="#B2182B", width=0.72, label="Up in DN")
axE.set_xticks(xx); axE.set_xticklabels([z.replace("Zone ", "") for z in order_z], fontsize=7)
axE.yaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))   # no out-of-view ticks
axE.set(xlabel="Zone", ylabel="# DEGs")
axE.legend(fontsize=7, loc="upper right", handlelength=0.9)

common.save_fig(fig, "Fig2_zonation")
print("saved Fig2_zonation")
