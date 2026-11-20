"""SFig 1: NMF metaprograms (183 mm). a MP heatmap, b selected MP boxplots, c pathway bubble."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import common
from common import SD, FIGDIR

common.set_nc_style()
mp = pd.read_csv(SD / "mp_scores.csv", comment="#").set_index("sample")
mp_de = pd.read_csv(SD / "mp_de.csv", comment="#")
pw = pd.read_csv(SD / "mp_pathways.csv", comment="#")

GRP_COL = {"NC": "#999999", "DN1": "#0072B2", "DN2": "#D55E00"}

fig = plt.figure(figsize=(common.W_DOUBLE, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.6, height_ratios=[1, 1.15])

# --- a: MP heatmap (z-scored per MP, samples grouped by group)
axA = fig.add_subplot(gs[0, 0]); common.panel_label(axA, "A")
mpcols = [c for c in mp.columns if c != "group1"]
order_s = mp.sort_values("group1").index.tolist()
mat = mp.loc[order_s, mpcols]
zmat = (mat - mat.mean(0)) / mat.std(0)
vmax = np.nanmax(np.abs(zmat.values)) or 1
im = axA.imshow(zmat.T.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
axA.set_xticks(range(len(order_s))); axA.set_xticklabels(order_s, fontsize=6.5, rotation=90)
axA.legend(handles=[Patch(fc=GRP_COL[g], label=g) for g in ["NC", "DN1", "DN2"]],
           fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.08), ncol=3,
           handlelength=0.9, columnspacing=0.6)
for i, s in enumerate(order_s):
    axA.axvspan(i - 0.5, i + 0.5, ymin=0.98, ymax=1.02, color=GRP_COL[mp.loc[s, "group1"]], clip_on=False)
axA.set_yticks(range(len(mpcols))); axA.set_yticklabels(mpcols, fontsize=7.5)
cb = plt.colorbar(im, ax=axA, shrink=0.7, pad=0.02); cb.set_label("z-scored MP score", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)

# --- b: selected MP boxplots with Wilcoxon p
axB = fig.add_subplot(gs[0, 1]); common.panel_label(axB, "B")
de_g2 = mp_de[mp_de.comparison == "Group2_DN-VS-NC"].set_index("mp")
sel = ["MP2", "MP6", "MP7", "MP1", "MP9"]
for i, m in enumerate(sel):
    xx = i * 3 + np.array([0, 1, 2])
    data = [mp[mp.group1 == g][m].values for g in ["NC", "DN1", "DN2"]]
    bp = axB.boxplot(data, positions=xx, widths=0.85, patch_artist=True,
                     medianprops=dict(color="black", lw=0.8), flierprops=dict(markersize=1.5))
    for patch, g in zip(bp["boxes"], ["NC", "DN1", "DN2"]):
        patch.set_facecolor(GRP_COL[g]); patch.set_alpha(0.8)
    rng = np.random.default_rng(42)
    for xi, g in enumerate(["NC", "DN1", "DN2"]):
        yy = mp[mp.group1 == g][m].values
        axB.scatter(rng.normal(xx[xi], 0.06, len(yy)), yy, s=4, c="#444444",
                    alpha=0.35, edgecolors="none", zorder=3)
    if m in de_g2.index:
        p = de_g2.loc[m, "pvalue"]
        axB.text(np.mean(xx), 0.98, f"p={p:.2g}", ha="center", va="top",
                 transform=axB.get_xaxis_transform(), fontsize=6)
axB.set_xticks([i * 3 + 1 for i in range(5)]); axB.set_xticklabels(sel, fontsize=7.5)
axB.set_ylabel("MP score"); axB.tick_params(labelsize=6.5)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(fc=GRP_COL[g], label=g) for g in ["NC", "DN1", "DN2"]],
           fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3,
           handlelength=0.9, columnspacing=0.6)

# --- c: pathway bubble (top 8 pathways for MP1/MP3/MP5/MP6/MP10)
axC = fig.add_subplot(gs[1, :]); common.panel_label(axC, "C")
sel_mp = ["MP1", "MP3", "MP5", "MP6", "MP10"]
dd = pw[pw.group.isin(sel_mp)].copy()
dd = dd.dropna(subset=["padj", "overlap"])
dd = dd.sort_values("padj").groupby("group").head(6)
paths = dd.pathway.str.replace("KEGG_", "").str.replace("_", " ").str[:11].unique().tolist()
dd["path"] = dd.pathway.str.replace("KEGG_", "").str.replace("_", " ").str[:11]
yy = {p: i for i, p in enumerate(paths)}
sc = axC.scatter([sel_mp.index(dd.group.iloc[i]) for i in range(len(dd))],
                 [yy[dd.path.iloc[i]] for i in range(len(dd))],
                 s=dd.overlap.values * 4, c=-np.log10(dd.padj.values), cmap="magma",
                 edgecolors="#333333", linewidths=0.2)
axC.set_xticks(range(len(sel_mp))); axC.set_xticklabels(sel_mp, fontsize=7.5)
axC.set_yticks(list(yy.values())); axC.set_yticklabels(list(yy.keys()), fontsize=6.5)
cb2 = plt.colorbar(sc, ax=axC, shrink=0.7, pad=0.02); cb2.set_label("-log10 padj", fontsize=7.5); cb2.ax.tick_params(labelsize=6.5)
# bubble size legend (overlap = number of overlapping genes)
from matplotlib.lines import Line2D
size_handles = [Line2D([0], [0], marker="o", linestyle="", markersize=np.sqrt(o * 4) * 0.6,
                      markerfacecolor="#888888", label=f"{o}") for o in [2, 5, 10]]
axC.legend(handles=size_handles, title="Overlap (genes)", title_fontsize=7,
           fontsize=6.5, loc="lower left", bbox_to_anchor=(1.02, 0.0),
           handletextpad=0.5)

common.save_fig(fig, "SFig1_NMF")
print("saved SFig1_NMF")
