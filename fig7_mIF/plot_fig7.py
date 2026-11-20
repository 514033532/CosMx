"""Fig 7 composer — multipanel-skill architecture: one independent subfigure per panel
with constrained_layout (no overlap by construction). A: local views (hand-exported
crops); B: glomerular density boxplots; C: whole-slide fold-change lollipops."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "fig7_panels"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common
import panel_A, panel_B, panel_C, panel_D

common.set_nc_style()

def add_letters(panels, size=14):
    for letter, sf in panels.items():
        sf.text(0.0, 1.0, letter, transform=sf.transSubfigure,
                fontsize=size, fontweight="bold", va="top", ha="left")

fig = plt.figure(figsize=(common.W_DOUBLE, 9.7), layout="constrained")
fig.get_layout_engine().set(rect=(0.012, 0, 0.988, 0.985))
sfs = fig.subfigures(3, 1, height_ratios=[2.65, 0.5, 0.5], hspace=0.08)

# --- A: local views (crops, 3 regions x 4 slides)
sfA = sfs[0]
axA = sfA.subplots()
panel_B.draw(axA)          # panel_B module draws the crop grid onto one axes

# --- B: glomerular density boxplots (five subplots)
sfB = sfs[1]
bB = sfB.subfigures(1, 5, wspace=0.10)
axes_b = [bB[i].subplots() for i in range(5)]
panel_C.draw(axes_b)

# --- C: whole-slide fold-change lollipops
sfC = sfs[2]
axC = sfC.subplots()
panel_D.draw(axC)

common.panel_label(axA, "A")
common.panel_label(axes_b[0], "B")
common.panel_label(axC, "C")

common.save_fig(fig, "Fig7_mIF_validation")
print("saved Fig7_mIF_validation")
