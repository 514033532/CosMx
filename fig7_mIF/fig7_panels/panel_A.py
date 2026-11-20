"""Fig 7 panel A: whole-slide mIF composites (4 slides, true-scale 5-mm bars)."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import MANI

DISEASE, CONTROL = "#D55E00", "#0072B2"
SLIDE_ORDER = [("13号片", "disease"), ("3号片", "disease"), ("对照1", "control"), ("对照2", "control")]
SLIDE_LAB = {"13号片": "DKD #1", "3号片": "DKD #2", "对照1": "Control #1", "对照2": "Control #2"}

def load_rgb(path, max_h=1400):
    from PIL import ImageEnhance
    im = Image.open(path).convert("RGB")
    orig_h = im.height
    if im.height > max_h:
        im = im.resize((int(im.width * max_h / im.height), max_h))
    im = ImageEnhance.Brightness(im).enhance(1.25)      # uniform, whole-image
    im = ImageEnhance.Contrast(im).enhance(1.15)        # uniform, whole-image
    return np.asarray(im), 4.0 * orig_h / im.height

def draw(axes):
    for ax, (slide, grp) in zip(axes, SLIDE_ORDER):
        row = MANI[(MANI.slide == slide) & (MANI.panel == "composite")]
        arr, um_per_disp_px = load_rgb(row.path.iloc[0])
        ax.imshow(arr)
        ax.set_title(SLIDE_LAB[slide], fontsize=8, pad=4,
                     color=DISEASE if grp == "disease" else CONTROL)
        ax.axis("off"); ax.set_anchor("N")
        W, H = arr.shape[1], arr.shape[0]
        len_px = 5000 / um_per_disp_px
        ax.plot([W - 40 - len_px, W - 40], [H * 0.04, H * 0.04],
                color="white", lw=2.5, solid_capstyle="butt")
        ax.text(W - 40 - len_px / 2, H * 0.04 + H * 0.012, "5 mm",
                color="white", fontsize=7, ha="center", va="bottom")

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 3.6))
    axes = [fig.add_axes(common.inset_rect(fig, [0.06 + i * 0.23, 0.10, 0.21, 0.80])) for i in range(4)]
    draw(axes)
    common.save_panel(fig, "fig7", "panel_A")
    print("saved panel_A")
