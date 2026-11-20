"""Compose Fig 6 from the finished standalone panel images (multipanel Option 2).
Faithful to the reviewed panel PNGs (native resolution, text unchanged); gaps 12 mm.
NOTE: the composite is ~390 mm tall — a fidelity copy for review. The submission
composite (<=247 mm) is Fig6_communication produced by fig6.py (Option 1, draw
functions shared with the panels)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread
import common

MM = 25.4
MARGIN, GAP = 10, 12
W = 180 - 2 * MARGIN

panels = {
    "A": "fig6/panel_A.png",
    "B": "fig6/panel_B.png",
    "C": "fig6/panel_C.png",
    "D": "fig6/panel_D.png",
    "E": "fig6/panel_E.png",
}
base = common.FIGDIR.parent / "figures_panels"
imgs = {k: imread(base / p) for k, p in panels.items()}

def h_mm(img, w_mm):
    return w_mm * img.shape[0] / img.shape[1]

# row heights at target widths
hA = h_mm(imgs["A"], W)
hB = h_mm(imgs["B"], W)
wCD = (W - GAP) / 2
hC = h_mm(imgs["C"], wCD)
hD = h_mm(imgs["D"], wCD)
hE = h_mm(imgs["E"], W)
H = MARGIN + hA + GAP + hB + GAP + max(hC, hD) + GAP + hE + MARGIN

fig = plt.figure(figsize=(180 / MM, H / MM), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")

def place(img, x_mm, y_top_mm, w_mm):
    h = h_mm(img, w_mm)
    ax.imshow(img, extent=[x_mm / 180, (x_mm + w_mm) / 180,
                           (H - y_top_mm - h) / H, (H - y_top_mm) / H])

y = MARGIN
place(imgs["A"], MARGIN, y, W); y += hA + GAP
place(imgs["B"], MARGIN, y, W); y += hB + GAP
hrow = max(hC, hD)
place(imgs["C"], MARGIN, y + (hrow - hC) / 2, wCD)
place(imgs["D"], MARGIN + wCD + GAP, y + (hrow - hD) / 2, wCD)
y += hrow + GAP
place(imgs["E"], MARGIN, y, W)

# panel letters at each block's top-left
for letter, (xx, yy) in {"A": (MARGIN, MARGIN), "B": (MARGIN, MARGIN + hA + GAP),
                          "C": (MARGIN, MARGIN + hA + GAP + hB + GAP),
                          "D": (MARGIN + wCD + GAP, MARGIN + hA + GAP + hB + GAP),
                          "E": (MARGIN, H - MARGIN - hE)}.items():
    ax.text(xx / 180, (H - yy) / H, letter, fontsize=14, fontweight="bold",
            va="top", ha="left", color="#333333")

out = common.FIGDIR / "Fig6_communication_images"
common.save_fig(fig, "Fig6_communication_images")
print(f"saved Fig6_communication_images ({H:.0f} mm tall — review copy)")
