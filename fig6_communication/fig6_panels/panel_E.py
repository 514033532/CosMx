"""Fig 6 panel E: key single ligand-receptor pairs (per-sample probabilities, NC vs DN)."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from data import LS, LV, CMP

KEY_PAIRS = [
 ("Distal tubular cell -> Neutrophil SPP1 -> ITGAV_ITGB5", "SPP1→ITGAV_ITGB5", "DT→Neutro", "up"),
 ("Cortical fibroblasts -> Distal tubular cell EGF -> EGFR", "EGF→EGFR", "cFibro→DT", "up"),
 ("Distal tubular cell -> Monocyte cell KLK3 -> NGFR", "KLK3→NGFR", "DT→Mono", "up"),
 ("Podocytes -> Podocytes VEGFA -> FLT1", "VEGFA→FLT1", "Podo→Podo", "down"),
 ("Podocytes -> Vascular Endothelial Cell VEGFA -> FLT1", "VEGFA→FLT1", "Podo→VEC", "down"),
 ("Podocytes -> M2 Macrophage PDGFA -> PDGFRB", "PDGFA→PDGFRB", "Podo→M2", "down")]

def draw(axes):
    rng = np.random.default_rng(42)
    for i, ((pair, lig_label, cell_label, direction), axe) in enumerate(zip(KEY_PAIRS, axes)):
        sub = LS[(LS.comparison == CMP) & (LS.lr_pair == pair)]
        data = [sub[sub.group1 == g].value.values for g in ["NC", "DN1", "DN2"]]
        dn_vals = np.concatenate([data[1], data[2]]) if len(data[1]) and len(data[2]) else data[1]
        data2 = [data[0], dn_vals]
        bp = axe.boxplot(data2, positions=[0, 1], widths=0.5, patch_artist=True,
                         medianprops=dict(color="black", lw=1), flierprops=dict(markersize=0, alpha=0))
        for patch, c in zip(bp["boxes"], ["#0072B2", "#D55E00"]):
            patch.set_facecolor(c); patch.set_alpha(0.75)
        for xi, d in enumerate(data2):
            if len(d):
                axe.scatter(rng.normal(xi, 0.07, len(d)), d, s=8, c="#444444",
                            alpha=0.6, edgecolors="none", zorder=3)
        axe.set_yscale("log")
        axe.yaxis.set_major_locator(plt.LogLocator(base=10, numticks=3))
        p = LV[(LV.comparison == CMP) & (LV.lr_pair == pair)].pvalue.values
        if len(p):
            ymax = max([max(d) for d in data2 if len(d)])
            axe.text(0.02, 1.06, f"p = {p[0]:.2g}", ha="left", va="bottom",
                     transform=axe.transAxes, fontsize=6.5,
                     bbox=dict(facecolor="white", edgecolor="none", pad=1))
        axe.set_xticks([0, 1]); axe.set_xticklabels(["NC", "DN"], fontsize=7)
        for t in list(axe.get_xticklabels()) + list(axe.get_yticklabels()):
            t.set_fontweight("bold")
        axe.set_xlabel(f"{lig_label}\n{cell_label}", fontsize=7, fontweight="bold", linespacing=1.2)
        if i == 0:
            axe.set_ylabel("Comm. prob. (log)", fontsize=8, fontweight="bold")
        else:
            axe.set_yticklabels([])
        axe.set_title("Up in DKD" if direction == "up" else "Down in DKD",
                      fontsize=7, color="#B2182B" if direction == "up" else "#2166AC",
                      loc="right", pad=2)

if __name__ == "__main__":
    common.set_nc_style()
    fig = plt.figure(figsize=(common.W_DOUBLE, 3.0))
    axes = [fig.add_axes(common.inset_rect(fig, [0.10 + i * 0.148, 0.18, 0.13, 0.62])) for i in range(6)]
    draw(axes)
    fig.text(0.5, 0.02, "Group", ha="center", fontsize=7.5, fontweight="bold")
    common.save_panel(fig, "fig6", "panel_E")
    print("saved panel_E")
