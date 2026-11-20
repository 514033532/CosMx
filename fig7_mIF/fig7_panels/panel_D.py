"""Fig 7 panel D: whole-slide density fold changes (lollipops, descriptive)."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import FC

DISEASE, CONTROL = "#D55E00", "#0072B2"
MARKERS = ["cd8", "cd20", "cd4", "cd68", "ctla4"]
MARKER_LAB = {"cd8": "CD8", "cd20": "CD20", "cd4": "CD4", "cd68": "CD68", "ctla4": "CTLA4"}

def draw(ax):
    fcd = FC[FC.metric == "dens"].set_index("marker").loc[MARKERS].sort_values("fold_change", ascending=False)
    yy = np.arange(len(fcd))[::-1]
    cols = [DISEASE if v > 1 else CONTROL for v in fcd.fold_change.values]
    ax.scatter(fcd.fold_change.values, yy, c=cols, s=34, zorder=3, edgecolors="#333333", linewidths=0.3)
    for yi, v in enumerate(fcd.fold_change.values):
        ax.plot([1, v], [yy[yi], yy[yi]], color="#999999", lw=1, zorder=1)
        ax.text(v, yy[yi] + 0.14, f"{v:.2f}×", fontsize=7, va="bottom", ha="center")
    ax.axvline(1, color="#333333", lw=0.8, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels([MARKER_LAB[m] for m in fcd.index], fontsize=7.5)
    ax.set_xlabel("Whole-slide density fold change (DN / control)")
    ax.set_xlim(0, 6.2)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc=DISEASE, label="Increased in DKD"),
                       Patch(fc=CONTROL, label="Decreased in DKD")],
              fontsize=7, loc="lower right", handlelength=0.9)
    ax.text(0.02, 0.97, "n = 2 slides/group (descriptive)", transform=ax.transAxes,
            fontsize=7, ha="left", va="top", color="#555555")

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 2.2))
    ax = fig.add_axes(common.inset_rect(fig, [0.14, 0.16, 0.82, 0.74]))
    draw(ax)
    common.save_panel(fig, "fig7", "panel_D")
    print("saved panel_D")
