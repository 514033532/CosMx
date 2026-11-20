"""Fig 5: three-dimensional differential expression (183 mm). Collision-free gene labels."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import common
from common import SD, FIGDIR

common.set_nc_style()
ct_de = pd.read_csv(SD / "ct_prop_de.csv", comment="#")
deg = pd.read_csv(SD / "ct_deg.csv", comment="#")
cent = pd.read_csv(SD / "centrality_de.csv", comment="#")

CMP = "Group2_DN-VS-NC"
FDR_T, FC_T = 0.05, 0.25
UP, DN, NS = "#B2182B", "#2166AC", "#BBBBBB"
DOT = dict(s=3, linewidths=0)
deg = deg[deg.comparison == CMP]

def dim_genes(celltypes):
    sub = deg[deg.celltype.isin(celltypes)]
    g = sub.groupby("gene").agg(mean_logFC=("logFC", "mean"), min_fdr=("fdr", "min")).reset_index()
    g["nlogp"] = -np.log10(g["min_fdr"].replace(0, 1e-300))
    g["sig"] = np.where((g.min_fdr < FDR_T) & (g.mean_logFC > FC_T), "up",
               np.where((g.min_fdr < FDR_T) & (g.mean_logFC < -FC_T), "down", "ns"))
    return g

prop_sig = set(ct_de[(ct_de.comparison == CMP) & (ct_de.significance == "Sig")].celltype)
ga = dim_genes(prop_sig)
gb = dim_genes(set(deg.celltype))
cent_d = cent[(cent.comparison == CMP) & (cent.metric == "Degree")].copy()
thr = cent_d.mean_diff.abs().quantile(0.75)
cent_sig = set(cent_d[cent_d.mean_diff.abs() >= thr].celltype)
gc = dim_genes(cent_sig)

fig = plt.figure(figsize=(common.W_DOUBLE, 3.8))
gs = fig.add_gridspec(1, 3, wspace=0.38)

def draw(ax, g, title, legend=False):
    ax.scatter(g[g.sig == "ns"].mean_logFC, g[g.sig == "ns"].nlogp, c=NS, rasterized=True, **DOT)
    ax.scatter(g[g.sig == "up"].mean_logFC, g[g.sig == "up"].nlogp, c=UP, rasterized=True, **DOT, label="Up in DN")
    ax.scatter(g[g.sig == "down"].mean_logFC, g[g.sig == "down"].nlogp, c=DN, rasterized=True, **DOT, label="Down in DN")
    top = g[g.sig != "ns"].reindex(g[g.sig != "ns"].mean_logFC.abs().sort_values(ascending=False).index).head(10)
    common.deoverlap_labels(ax, top.mean_logFC.values, top.nlogp.values, top.gene.values, fontsize=7)
    ax.set_title(title, fontsize=8)
    ax.set(xlabel="Mean log2 FC", ylabel="−log10 min FDR")
    ax.xaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))   # sparse, in-range ticks
    for xv in (FC_T, -FC_T):
        ax.axvline(xv, color="#999999", lw=0.6)
    ax.axhline(-np.log10(FDR_T), color="#999999", lw=0.6)
    if legend:
        ax.legend(fontsize=7, loc="upper right", handlelength=0.9, handletextpad=0.3)

titles = ["Proportion dimension", "Cluster dimension", "Centrality dimension"]
for i, (g, t) in enumerate([(ga, titles[0]), (gb, titles[1]), (gc, titles[2])]):
    ax = fig.add_subplot(gs[0, i]); common.panel_label(ax, chr(97 + i))
    draw(ax, g, t, legend=(i == 0))

common.save_fig(fig, "Fig5_3D_DE")
print("saved Fig5_3D_DE")
