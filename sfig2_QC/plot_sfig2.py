"""SFig 2: QC & robustness (183 mm). a QC overview, b threshold sensitivity, c coexpression overlap."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import common
from common import SD, FIGDIR

common.set_nc_style()
qc = pd.read_csv(SD / "qc_overview.csv", comment="#")
sens = pd.read_csv(SD / "mif_sensitivity.csv", comment="#")
co = pd.read_csv(SD / "mif_coexpression.csv", comment="#")

fig = plt.figure(figsize=(common.W_DOUBLE, 7.2))
fig.subplots_adjust(left=0.16)   # room for QC metric labels in panel a
gs = fig.add_gridspec(3, 1, hspace=0.42, height_ratios=[1.1, 1.0, 1.2])

# --- a: QC overview bars
axA = fig.add_subplot(gs[0, 0]); common.panel_label(axA, "A")
qmap = {"Total_Cells_Num": "Cells", "Total_FOVs_Num": "FOVs", "Mean_Cell_Area": "Cell area (µm²)",
        "Mean_Exp_per_Cell": "Transcripts / cell", "Mean_UniqueGene_per_Cell": "Unique genes / cell",
        "Mean_NegativeProbe": "Neg. probes / cell", "Mean_FalseCode": "False codes / cell"}
qq = qc[qc.metric.isin(qmap)].copy()
qq["label"] = qq.metric.map(qmap); qq["value"] = pd.to_numeric(qq.value, errors="coerce")
qq = qq.sort_values("value", ascending=False)
axA.barh(np.arange(len(qq))[::-1], qq.value.values, color="#4E79A7", height=0.6)
axA.set_yticks(np.arange(len(qq))[::-1]); axA.set_yticklabels(qq.label.values, fontsize=6.5)
for yi, v in enumerate(qq.value.values):
    axA.text(v, len(qq) - 1 - yi, f" {v:g}", fontsize=6.5, va="center")
axA.set_xlabel("value (log10 scale)"); axA.set_xscale("log")

# --- b: threshold sensitivity lollipops
axB = fig.add_subplot(gs[1, 0]); common.panel_label(axB, "B")
markers = sens.marker.unique().tolist()
thrs = sorted(sens.threshold.unique())
for i, m in enumerate(markers):
    sub = sens[sens.marker == m].sort_values("threshold")
    xx = i * 3 + np.arange(len(sub))
    cols = ["#B2182B" if v > 1 else "#2166AC" for v in sub.fold_change]
    axB.scatter(xx, sub.fold_change.values, c=cols, s=14, zorder=3, edgecolors="#333333", linewidths=0.2)
    for xi, (t, fc) in enumerate(zip(sub.threshold, sub.fold_change)):
        axB.plot([i * 3 + xi, i * 3 + xi], [1, fc], color="#999999", lw=0.6, zorder=1)
axB.axhline(1, color="#333333", lw=0.8, ls="--")
axB.set_xticks([i * 3 + 1 for i in range(len(markers))])
axB.set_xticklabels([m.upper() for m in markers], fontsize=7.5, rotation=45)
axB.set_ylabel("Fold change (DN / control)")
axB.set_title("Otsu threshold ×0.8 / ×1.0 / ×1.2 (3 points per marker)", fontsize=8)
axB.tick_params(labelsize=6.5)
# red/blue legend
from matplotlib.patches import Patch
axB.legend(handles=[Patch(fc="#B2182B", label="fold change > 1"),
                    Patch(fc="#2166AC", label="fold change < 1")],
           fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5), handlelength=0.9)

# --- c: coexpression overlap heatmap
axC = fig.add_subplot(gs[2, 0]); common.panel_label(axC, "C")
co["pair"] = co.m1.str.upper() + "+" + co.m2.str.upper()
mat = co.pivot_table(index="pair", columns="image_id", values="pct_of_m1", aggfunc="first")
mat = mat[["24R2590SLZA", "24R2889SLZA", "24R8601SLZA", "24R8119SLZA"]]
vmax = mat.values.max()
im = axC.imshow(mat.values, cmap="YlOrRd", vmin=0, vmax=vmax, aspect="auto")
axC.axvline(1.5, color="white", lw=2)   # separate DN from Ctrl columns
axC.set_xticks(range(4)); axC.set_xticklabels(["DN #1", "DN #2", "Ctrl #1", "Ctrl #2"], fontsize=7, rotation=30)
axC.set_yticks(range(len(mat))); axC.set_yticklabels(mat.index, fontsize=6.5)
for i in range(mat.shape[0]):
    for j in range(4):
        axC.text(j, i, f"{mat.values[i, j]:.0f}", ha="center", va="center", fontsize=6,
                 color="white" if mat.values[i, j] > vmax * 0.5 else "black")
cb = plt.colorbar(im, ax=axC, shrink=0.8, pad=0.02); cb.set_label("% of m1+ cells", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)

common.save_fig(fig, "SFig2_QC_robustness")
print("saved SFig2_QC_robustness")
