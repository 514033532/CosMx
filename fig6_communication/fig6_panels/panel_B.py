"""Fig 6 panel B: cell-pair LR-count DE heatmap, split into two stacked sub-panels
(8 + 7 pairs) so every column label is fully readable (~7.4 mm pitch)."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import LN

COMPS = ["Group1_DN1-VS-NC", "Group1_DN2-VS-NC", "Group2_DN-VS-NC"]
COMP_LAB = {"Group1_DN1-VS-NC": "DN1 vs NC", "Group1_DN2-VS-NC": "DN2 vs NC", "Group2_DN-VS-NC": "DN vs NC"}

def _sub(ax, mat, smat, vmax, labels):
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
    ax.set_yticks(range(3)); ax.set_yticklabels(mat.index, fontsize=8, fontweight="bold")
    for i in range(3):
        for j in range(len(labels)):
            v = mat.values[i, j]
            if pd.notna(v):
                ax.text(j, i, "*" if smat.values[i, j] == "Sig" else "", ha="center", va="center",
                        fontsize=7, color="white" if abs(v) > vmax * 0.55 else "black")
    return im

def draw(ax_top, ax_bot):
    """Top: pairs 1-8; bottom: pairs 9-15. Shared colour scale."""
    bb = LN[LN.comparison.isin(COMPS)].copy()
    bb["cmp"] = bb.comparison.map(COMP_LAB)
    top = bb.groupby("cell_pair")["mean_diff"].agg(lambda s: s.abs().max()).sort_values(ascending=False).head(15).index.tolist()
    parts = [top[:8], top[8:]]
    vmax = np.nanmax(np.abs(bb[bb.cell_pair.isin(top)].pivot_table(
        index="cmp", columns="cell_pair", values="mean_diff", aggfunc="first").values)) or 1
    ims = []
    for ax, pairs in zip((ax_top, ax_bot), parts):
        mat = bb.pivot_table(index="cmp", columns="cell_pair", values="mean_diff",
                             aggfunc="first").reindex(index=list(COMP_LAB.values()), columns=pairs)
        smat = bb.pivot_table(index="cmp", columns="cell_pair", values="significance",
                              aggfunc="first").reindex(index=list(COMP_LAB.values()), columns=pairs)
        ims.append(_sub(ax, mat, smat, vmax, [common.pair_label(q) for q in pairs]))
    return ims[0]

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 5.2))
    ax_top = fig.add_axes(common.inset_rect(fig, [0.14, 0.62, 0.70, 0.30]))
    ax_bot = fig.add_axes(common.inset_rect(fig, [0.14, 0.10, 0.70, 0.42]))
    im = draw(ax_top, ax_bot)
    cb = plt.colorbar(im, ax=ax_bot, fraction=0.05, pad=0.02)
    cb.set_label("Δ LR count", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    common.save_panel(fig, "fig6", "panel_B")
    print("saved panel_B")
