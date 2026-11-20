"""Fig 7 panel B: local views — hand-exported crops, 3 regions x 4 slides."""
from pathlib import Path
import os
import sys
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common

DISEASE, CONTROL = "#D55E00", "#0072B2"
SLIDE_ORDER = [("13号片", "disease"), ("3号片", "disease"), ("对照1", "control"), ("对照2", "control")]
SLIDE_LAB = {"13号片": "DKD #1", "3号片": "DKD #2", "对照1": "Control #1", "对照2": "Control #2"}
CROPS = os.environ.get(
    "MIF_CROPS_DIR",
    str(Path(__file__).resolve().parents[4] / "panel2_analysis" / "result20261020"),
)
SLIDE_NUM = {"13号片": "2590", "3号片": "2889", "对照1": "8601", "对照2": "8119"}
# distinct ROI per slide so adjacent control columns do not look duplicated
CROP_IDX = {"13号片": "-1", "3号片": "-2", "对照1": "-4", "对照2": "-5"}
REGIONS = [("肾小球（Glomerulus）", "Glomerulus"),
           ("肾小管间质（Interstitium）", "Interstitium"),
           ("perivascular", "Perivascular")]

def draw(ax):
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    n_reg, n_sl = len(REGIONS), len(SLIDE_ORDER)
    left_margin, top_margin = 0.075, 0.06
    gap_x, gap_y = 0.035, 0.08   # >= 6 mm at 180 mm width / 80 mm height
    cell_w = (1 - left_margin - 0.005 - gap_x * (n_sl - 1)) / n_sl
    cell_h = (1 - top_margin - 0.005 - gap_y * (n_reg - 1)) / n_reg
    for j, (slide, grp) in enumerate(SLIDE_ORDER):
        xc = left_margin + j * (cell_w + gap_x) + cell_w / 2
        ax.text(xc, 0.995, SLIDE_LAB[slide], ha="center", va="top", fontsize=7.5,
                color=DISEASE if grp == "disease" else CONTROL, fontweight="bold")
    ax.text(0.5, 0.0, "Local views: merge of DAPI / CD20 / CD68 / CD8 / CD4 (all channels)",
            ha="center", va="bottom", fontsize=7, color="#555555")
    for i, (cn, en) in enumerate(REGIONS):
        y_top = 1 - top_margin - i * (cell_h + gap_y)
        ax.text(0.0, y_top - cell_h / 2, en, ha="right", va="center", fontsize=7.5,
                rotation=90, color="#333333")
        for j, (slide, grp) in enumerate(SLIDE_ORDER):
            x0 = left_margin + j * (cell_w + gap_x)
            y0 = y_top - cell_h
            cand = glob.glob(f"{CROPS}/{cn}/{SLIDE_NUM[slide]}{CROP_IDX[slide]}.PNG")
        if not cand:
            cand = glob.glob(f"{CROPS}/{cn}/{SLIDE_NUM[slide]}-*.PNG")
            if cand:
                img = Image.open(sorted(cand)[0]).convert("RGB")
                # fit the image into the cell preserving its aspect ratio
                ih, iw = img.size[1], img.size[0]
                ar = iw / ih
                if cell_w / cell_h > ar:
                    w = cell_h * ar; x0c = x0 + (cell_w - w) / 2
                    ax.imshow(img, aspect="equal", extent=[x0c, x0c + w, y0, y_top])
                else:
                    h = cell_w / ar; y0c = y0 + (cell_h - h) / 2
                    ax.imshow(img, aspect="equal", extent=[x0, x0 + cell_w, y0c, y0c + h])

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 4.6))
    ax = fig.add_axes(common.inset_rect(fig, [0.02, 0.02, 0.96, 0.96]))
    draw(ax)
    common.save_panel(fig, "fig7", "panel_B")
    print("saved panel_B")
