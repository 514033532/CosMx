"""Fig 6 panel C: top up/down ligand-receptor pairs (16 bars, abbreviated labels)."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import LV, CMP

def parse_lr(s):
    """'SRC -> TGT LIG -> REC' -> (src, tgt, lig, rec)."""
    parts = s.split(" -> ")
    if len(parts) != 3:
        return None
    src, tgt_lig, rec = parts
    tgt, lig = None, None
    for ct in common.CELLTYPES:
        if tgt_lig.startswith(ct + " "):
            tgt, lig = ct, tgt_lig[len(ct) + 1:]
            break
    return (src, tgt, lig, rec) if tgt else None

def short_lr_label(s):
    parsed = parse_lr(s)
    if not parsed:
        return s[:58]
    src, tgt, lig, rec = parsed
    cell = f"{common.CT_ABBR[src]}→{common.CT_ABBR[tgt]}" if src != tgt else common.CT_ABBR[src]
    return f"{cell} {lig[:8]}→{rec[:8]}"

def draw(ax):
    vv = LV[LV.comparison == CMP].copy()
    vv = vv[vv.significance == "Sig"]
    up = vv[vv.change.str.contains("up")].sort_values("mean_diff", ascending=False).head(8)
    dn = vv[vv.change.str.contains("down")].sort_values("mean_diff").head(8)
    both = pd.concat([up, dn])
    yy = np.arange(len(both))[::-1]
    cols = ["#B2182B"] * len(up) + ["#2166AC"] * len(dn)
    ax.barh(yy, both["mean_diff"].values, color=cols, height=0.6)
    ax.set_yticks(yy)
    ax.set_yticklabels([short_lr_label(r.lr_pair) for r in both.itertuples()], fontsize=7)
    for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        t.set_fontweight("bold")
    ax.set_xlabel("Δ communication probability", fontsize=8, fontweight="bold")
    ax.xaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))
    ax.axvline(0, color="#333333", lw=1.0)
    ax.set_axisbelow(True)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc="#B2182B", label="Up in DKD"),
                       Patch(fc="#2166AC", label="Down in DKD")],
              fontsize=7, loc="lower right", handlelength=0.9)
    ax.set_title("Ligand–receptor pairs (DN vs NC, significant)", fontsize=7.5)

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 4.0))
    ax = fig.add_axes(common.inset_rect(fig, [0.22, 0.12, 0.74, 0.82]))
    draw(ax)
    common.save_panel(fig, "fig6", "panel_C")
    print("saved panel_C")
