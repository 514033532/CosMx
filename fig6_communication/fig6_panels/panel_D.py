"""Fig 6 panel D: communication centrality DE heatmap (degree/closeness/clustering)."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import CD, CMP

def draw(ax):
    dd = CD[CD.comparison == CMP].copy()
    mats = dd.pivot_table(index="metric", columns="celltype", values="mean_diff", aggfunc="first").reindex(columns=common.CELLTYPES, index=["Degree", "Closeness", "Clustering"])
    smats = dd.pivot_table(index="metric", columns="celltype", values="significance", aggfunc="first").reindex(columns=common.CELLTYPES, index=["Degree", "Closeness", "Clustering"])
    vmax = np.nanmax(np.abs(mats.values)) or 1
    im = ax.imshow(mats.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    CT_LAB = [common.CT_SHORT[c] for c in common.CELLTYPES]
    ax.set_xticks(range(18)); ax.set_xticklabels(CT_LAB, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(3)); ax.set_yticklabels(mats.index, fontsize=8)
    for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        t.set_fontweight("bold")
    for i in range(3):
        for j in range(18):
            v = mats.values[i, j]
            if pd.notna(v):
                ax.text(j, i, "*" if smats.values[i, j] == "Sig" else "", ha="center", va="center",
                        fontsize=7, color="white" if abs(v) > vmax * 0.55 else "black")
    ax.set_xlabel("Cell type", fontsize=8, fontweight="bold")
    return im

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 2.6))
    ax = fig.add_axes(common.inset_rect(fig, [0.14, 0.30, 0.78, 0.55]))
    im = draw(ax)
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Δ centrality", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    common.save_panel(fig, "fig6", "panel_D")
    print("saved panel_D")
