"""Fig 7 panel C: glomerular marker density boxplots (BH-adjusted p annotated)."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import GDENS, MWU

DISEASE, CONTROL = "#D55E00", "#0072B2"
MARKERS = ["cd8", "cd20", "cd4", "cd68", "ctla4"]
MARKER_LAB = {"cd8": "CD8", "cd20": "CD20", "cd4": "CD4", "cd68": "CD68", "ctla4": "CTLA4"}

def draw(axes):
    np.random.seed(42)
    mwu_g = MWU[MWU.roi_class == "Glomerulus"].set_index("marker")
    for i, (m, ax) in enumerate(zip(MARKERS, axes)):
        data = [GDENS[(GDENS.marker == m) & (GDENS.group == g)].density_mm2.values for g in ["control", "disease"]]
        bp = ax.boxplot(data, positions=[0, 1], widths=0.5, patch_artist=True,
                        medianprops=dict(color="black", lw=1), flierprops=dict(markersize=1.5, alpha=0.4))
        for patch, g, c in zip(bp["boxes"], ["control", "disease"], [CONTROL, DISEASE]):
            patch.set_facecolor(c); patch.set_alpha(0.85)
        for xi, g in enumerate(["control", "disease"]):
            yy = GDENS[(GDENS.marker == m) & (GDENS.group == g)].density_mm2.values
            ax.scatter(np.random.normal(xi, 0.06, len(yy)), yy, s=2.5, c="#444444", alpha=0.25, zorder=3, edgecolors="none")
        ax.set_title(MARKER_LAB[m], fontsize=8)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Ctrl", "DN"], fontsize=7.5)
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax.tick_params(labelsize=7)
        if m in mwu_g.index:
            p = mwu_g.loc[m, "p_adj_bh"]
            ax.text(0.02, 0.97, f"{common.sig_stars(p)} p={p:.2g}",
                    ha="left", va="top", transform=ax.transAxes, fontsize=7)

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 2.8))
    axes = [fig.add_axes(common.inset_rect(fig, [0.08 + i * 0.182, 0.18, 0.155, 0.68])) for i in range(5)]
    draw(axes)
    fig.text(0.5, 0.02, "Glomerular marker density (cells/mm²) — n = 600 ROIs/group, two-sided MWU, BH-adjusted p",
             ha="center", fontsize=7.5, color="#555555")
    common.save_panel(fig, "fig7", "panel_C")
    print("saved panel_C")
