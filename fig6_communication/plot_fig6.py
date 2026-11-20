"""Fig 6 composer — multipanel-skill architecture: one independent subfigure per panel
with constrained_layout, so every panel packs to its OWN labels and nothing can overlap
by construction. Panel letters sit in each subfigure's frame (skill helper pattern)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "fig6_panels"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common
import panel_A, panel_B, panel_C, panel_D, panel_E

common.set_nc_style()

def add_letters(panels, size=14):
    """Bold letter at each subfigure's own top-left (multipanel skill helper)."""
    for letter, sf in panels.items():
        sf.text(0.0, 1.0, letter, transform=sf.transSubfigure,
                fontsize=size, fontweight="bold", va="top", ha="left")

fig = plt.figure(figsize=(common.W_DOUBLE, 9.45), layout="constrained")
fig.get_layout_engine().set(rect=(0.012, 0, 0.988, 0.985))   # hair of margin for letters
sfs = fig.subfigures(4, 1, height_ratios=[1.15, 0.75, 1.45, 0.65], hspace=0.06)

# --- A: paired communication-strength matrices
sfA = sfs[0]
aA = sfA.subfigures(1, 2, wspace=0.10)
axA0 = aA[0].subplots(); axA1 = aA[1].subplots()
im_a = panel_A.draw(axA0, axA1)
cba = fig.colorbar(im_a, ax=axA1, fraction=0.03, pad=0.01)
cba.set_label("log1p strength", fontsize=7.5); cba.ax.tick_params(labelsize=7)

# --- B: side-by-side sub-heatmaps (pairs 1-8 | 9-15)
sfB = sfs[1]
bB = sfB.subfigures(1, 2, wspace=0.12, width_ratios=[1, 1])
axB_l = bB[0].subplots(); axB_r = bB[1].subplots()
im_b = panel_B.draw(axB_l, axB_r)
cbb = fig.colorbar(im_b, ax=axB_r, fraction=0.04, pad=0.01)
cbb.set_label("Δ LR count", fontsize=7.5); cbb.ax.tick_params(labelsize=7)

# --- C (LR bars) | D (centrality) side by side
sfCD = sfs[2]
cD = sfCD.subfigures(1, 2, wspace=0.12, width_ratios=[1, 1.25])
axC = cD[0].subplots()
panel_C.draw(axC)
axD = cD[1].subplots()
im_d = panel_D.draw(axD)
cbd = fig.colorbar(im_d, ax=axD, fraction=0.04, pad=0.01)
cbd.set_label("Δ centrality", fontsize=7.5); cbd.ax.tick_params(labelsize=7)

# --- E: six mini boxplots
sfE = sfs[3]
eE = sfE.subfigures(1, 6, wspace=0.10)
axes_e = [eE[i].subplots() for i in range(6)]
panel_E.draw(axes_e)

common.panel_label(axA0, "A")
common.panel_label(axB_l, "B")
common.panel_label(axC, "C")
common.panel_label(axD, "D")
common.panel_label(axes_e[0], "E")

common.save_fig(fig, "Fig6_communication")
print("saved Fig6_communication")
