"""Fig 3: molecular pathology of zones 0/3/6 (183 mm). Collision-free gene labels,
solid hairline thresholds, clean pathway bars."""
import pandas as pd
import numpy as np
import textwrap
import matplotlib.pyplot as plt
import common
from common import SD, FIGDIR

common.set_nc_style()
zd = pd.read_csv(SD / "zone_deg.csv", comment="#")
zp = pd.read_csv(SD / "zone_pathways.csv", comment="#")
zd = zd[zd.comparison == "Group2_DN-VS-NC"]
zp = zp[zp.comparison == "Group2_DN-VS-NC"]

FDR_T, FC_T = 0.05, 0.25
ZONES = ["Zone 0", "Zone 3", "Zone 6"]
UP, DN, NS = "#B2182B", "#2166AC", "#BBBBBB"
DOT = dict(s=3, linewidths=0)

def wrap2(d):
    """KEGG description on <=2 lines of <=13 chars, ellipsis if longer (keeps labels
    inside their own column — never crossing into the neighbouring panel)."""
    lines = textwrap.wrap(d, 13)
    if len(lines) > 2:
        lines = lines[:2]; lines[1] += "…"
    return "\n".join(lines)

fig = plt.figure(figsize=(common.W_DOUBLE, 7.6))
fig.subplots_adjust(left=0.17)   # room for wrapped pathway names
gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.55, width_ratios=[1, 1, 1.12])

def volcano(ax, zdf, zone, legend=False):
    dd = zdf[zdf.zone == zone].copy()
    dd["sig"] = np.where((dd.fdr < FDR_T) & (dd.logFC > FC_T), "up",
                np.where((dd.fdr < FDR_T) & (dd.logFC < -FC_T), "down", "ns"))
    dd["nlogp"] = -np.log10(dd.fdr.replace(0, 1e-300))
    ax.scatter(dd[dd.sig == "ns"].logFC, dd[dd.sig == "ns"].nlogp, c=NS, rasterized=True, **DOT)
    ax.scatter(dd[dd.sig == "up"].logFC, dd[dd.sig == "up"].nlogp, c=UP, rasterized=True, **DOT, label="Up in DN")
    ax.scatter(dd[dd.sig == "down"].logFC, dd[dd.sig == "down"].nlogp, c=DN, rasterized=True, **DOT, label="Down in DN")
    top = dd[dd.sig != "ns"].reindex(dd[dd.sig != "ns"].logFC.abs().sort_values(ascending=False).index).head(8)
    common.deoverlap_labels(ax, top.logFC.values, top.nlogp.values, top.gene.values, fontsize=7)
    ax.set_title(zone, fontsize=8)
    ax.set(xlabel="log2 fold change", ylabel="−log10 FDR")
    ax.xaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))   # sparse, in-range ticks
    for xv in (FC_T, -FC_T):
        ax.axvline(xv, color="#999999", lw=0.6)
    ax.axhline(-np.log10(FDR_T), color="#999999", lw=0.6)
    if legend:
        ax.legend(fontsize=7, loc="upper right", handlelength=0.9, handletextpad=0.3)

def pathways(ax, pdf, zone):
    dd = pdf[pdf.zone == zone].copy()
    dd = dd.dropna(subset=["p.adjust"]).sort_values("p.adjust").head(8)
    dd["nlogp"] = -np.log10(dd["p.adjust"].clip(lower=1e-300))
    cols = [UP if n > 0 else DN for n in dd["NES"]]
    ax.barh(np.arange(len(dd))[::-1], dd["nlogp"], color=cols, height=0.58)
    ax.set_yticks(np.arange(len(dd))[::-1])
    ax.set_yticklabels([wrap2(d) for d in dd["Description"]], fontsize=6.5, linespacing=1.1)
    ax.set_xlabel("−log10 adj. p")
    ax.set_axisbelow(True)

for i, z in enumerate(ZONES):
    ax = fig.add_subplot(gs[0, i]); common.panel_label(ax, chr(97 + i))
    volcano(ax, zd, z, legend=(i == 0))
for i, z in enumerate(ZONES):
    ax = fig.add_subplot(gs[1, i]); common.panel_label(ax, chr(100 + i))
    pathways(ax, zp, z)
    if i == 2:
        ax.set_title("KEGG pathways (top 8 by adj. p)", fontsize=8)   # red/blue direction: panel a legend

common.save_fig(fig, "Fig3_zone_pathology")
print("saved Fig3_zone_pathology")
