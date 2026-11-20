"""Fig 6 panel A: mean communication strength matrices NC vs DN (18 x 18, log1p)."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import CS

def draw(ax0, ax1):
    ct_order = common.CELLTYPES
    nc = CS[CS.group1 == "NC"].groupby(["source", "target"])["strength"].mean().unstack("target").reindex(index=ct_order, columns=ct_order)
    dn = CS[CS.group1.isin(["DN1", "DN2"])].groupby(["source", "target"])["strength"].mean().unstack("target").reindex(index=ct_order, columns=ct_order)
    vmax = max(nc.values.max(), dn.values.max())
    im = ax0.imshow(np.log1p(nc.values), cmap="viridis", vmin=0, vmax=np.log1p(vmax), aspect="auto")
    ax1.imshow(np.log1p(dn.values), cmap="viridis", vmin=0, vmax=np.log1p(vmax), aspect="auto")
    CT_LAB = ["Podo","PEC","Vasc.endo","VSMC","PT","DT","Cort.fibro","Inf.fibro",
              "Act.fibro","Pericyte","M2 macro","Monocyte","Neutro","Mast","T","B","NK","Plasma"]
    for axx in (ax0, ax1):
        axx.set_xticks(range(18)); axx.set_xticklabels(CT_LAB, rotation=90, fontsize=7)
        axx.set_yticks(range(18)); axx.set_yticklabels(CT_LAB, fontsize=7)
        for t in list(axx.get_xticklabels()) + list(axx.get_yticklabels()):
            t.set_fontweight("bold")
        axx.text(0.5, -0.36, "Target", transform=axx.transAxes, ha="center",
                 va="top", fontsize=8, fontweight="bold")
    ax0.set_title("NC (n = 11)", fontsize=7.5)
    ax1.set_title("DN (n = 5)", fontsize=7.5)
    ax1.set_yticklabels([])
    ax1.tick_params(axis="y", length=0)          # no orphan tick marks on the right heatmap
    ax0.set_ylabel("Source", fontsize=8, fontweight="bold")
    return im

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 4.2))
    ax0 = fig.add_axes(common.inset_rect(fig, [0.16, 0.30, 0.35, 0.56]))
    ax1 = fig.add_axes(common.inset_rect(fig, [0.58, 0.30, 0.35, 0.56]))
    im = draw(ax0, ax1)
    cb = plt.colorbar(im, ax=ax1, fraction=0.05, pad=0.02)
    cb.set_label("log1p strength", fontsize=7.5); cb.ax.tick_params(labelsize=7)
    common.save_panel(fig, "fig6", "panel_A")
    print("saved panel_A")
