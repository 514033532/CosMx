"""NC manuscript style module — paths, palettes, helpers. All figure scripts import from here."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT  = Path("/home/zhh2/YUBO-202508/ranse/panel2_analysis/NC_manuscript_v1_20261009")
SD, FIGDIR, TABDIR = OUT/"source_data", OUT/"figures", OUT/"tables"
SRC  = Path("/home/zhh2/YUBO-202508/20251114_share")
MIF  = Path("/home/zhh2/YUBO-202508/ranse/panel2_analysis")
DRAFT = Path("/home/zhh2/YUBO-202508/jieguo-v1.txt")

W_SINGLE, W_DOUBLE = 88/25.4, 180/25.4   # inches — journal spec: 88 mm single, 180 mm double

# 18 main cell types — validator-optimized set (lightness band + chroma floor PASS;
# CVD worst adjacent pair 12.4 ΔE PASS after 2026-10-10 clarity pass (Neutrophil→dark green
# #1B7837, breaking the old #E0559E↔#C25B74 6.3-ΔE pair); normal-vision floor cannot reach
# ΔE≥15 for all 18 adjacent pairs — inherent to 18-color categorical encoding, mitigated by
# labeled legends everywhere + dotplot + per-type small-multiple UMAPs + group-level views.
# Full validator output in palette_check.txt.
CELLTYPE_COLORS = {
 "Podocytes": "#e6194B", "Parietal epithelial cell": "#f58231",
 "Vascular Endothelial Cell": "#4363d8", "VSMC": "#3cb44b",
 "Proximal tubular cell": "#00A196", "Distal tubular cell": "#9A6324",
 "Cortical fibroblasts": "#911eb4", "Inflammatory fibroblasts": "#f032e6",
 "Activated myofibroblasts": "#808000", "Pericyte": "#A1548C",
 "M2 Macrophage": "#9E2A2B", "Monocyte cell": "#7FA800",
 "Neutrophil": "#1B7837", "Mast": "#E0559E",
 "T": "#D4AF37", "B": "#00B3D6", "NK": "#5B3A8E", "Plasma": "#BC80BD",
}
MAJOR_GROUP = {
 "Podocytes": "Glomerular", "Parietal epithelial cell": "Glomerular",
 "Cortical fibroblasts": "Interstitial", "Inflammatory fibroblasts": "Interstitial",
 "Activated myofibroblasts": "Interstitial", "Pericyte": "Interstitial",
 "T": "Lymphoid", "B": "Lymphoid", "NK": "Lymphoid",
 "M2 Macrophage": "Macrophage",
 "Monocyte cell": "Myeloid", "Neutrophil": "Myeloid", "Mast": "Myeloid",
 "Plasma": "Plasma",
 "Proximal tubular cell": "Tubular", "Distal tubular cell": "Tubular",
 "Vascular Endothelial Cell": "Vascular", "VSMC": "Vascular",
}
GROUP_COLORS = {"Glomerular": "#e6194B", "Interstitial": "#911eb4",
 "Lymphoid": "#D4AF37", "Macrophage": "#9E2A2B", "Myeloid": "#7FA800",
 "Plasma": "#BC80BD", "Tubular": "#00A196", "Vascular": "#4363d8"}
ZONE_COLORS = {f"Zone {i}": c for i, c in enumerate(
 ["#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948",
  "#B07AA1","#FF9DA7","#9C755F","#BAB0AC","#86BCB6"])}
MIF_GROUP_COLORS = {"disease": "#D55E00", "control": "#0072B2"}
CELLTYPES = list(CELLTYPE_COLORS)
ZONES = [f"Zone {i}" for i in range(11)]

# thin dark edge for scatter markers so light colors stay visible on white
SCATTER_EDGE = {"edgecolors": "#333333", "linewidths": 0.1}

# short display names for 18 cell types (heatmap axes)
CT_SHORT = {"Podocytes": "Podocytes", "Parietal epithelial cell": "Parietal epithelial",
    "Vascular Endothelial Cell": "Vascular endothelial", "VSMC": "VSMC",
    "Proximal tubular cell": "Proximal tubular", "Distal tubular cell": "Distal tubular",
    "Cortical fibroblasts": "Cortical fibro.", "Inflammatory fibroblasts": "Inflammatory fibro.",
    "Activated myofibroblasts": "Act. myofibro.", "Pericyte": "Pericyte",
    "M2 Macrophage": "M2 macrophage", "Monocyte cell": "Monocyte", "Neutrophil": "Neutrophil",
    "Mast": "Mast", "T": "T", "B": "B", "NK": "NK", "Plasma": "Plasma"}

# compact nephrology-standard abbreviations (pair labels: PT/DT/PEC are standard usage)
CT_ABBR = {"Podocytes": "Podo", "Parietal epithelial cell": "PEC",
    "Vascular Endothelial Cell": "Vasc. endo.", "VSMC": "VSMC",
    "Proximal tubular cell": "PT", "Distal tubular cell": "DT",
    "Cortical fibroblasts": "cFibro", "Inflammatory fibroblasts": "iFibro",
    "Activated myofibroblasts": "aFibro", "Pericyte": "Pericyte",
    "M2 Macrophage": "M2 macro", "Monocyte cell": "Monocyte",
    "Neutrophil": "Neutro", "Mast": "Mast", "T": "T", "B": "B", "NK": "NK", "Plasma": "Plasma"}

MIN_FS = 7   # minimum font size anywhere in figures

def set_nc_style(font="Arial"):   # Arial probed available at /usr/share/fonts/truetype/msttcorefonts
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [font, "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "axes.linewidth": 1.0, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5,
        "lines.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
        "grid.color": "#E3E3E1", "grid.linewidth": 0.4,
        "figure.dpi": 300, "savefig.dpi": 300, "svg.fonttype": "none",
        "pdf.fonttype": 42, "legend.frameon": False, "legend.fontsize": 7,
        "axes.titlelocation": "left", "axes.titlesize": 7.5, "axes.titleweight": "normal"})

def inset_rect(fig, rect, margin_mm=3.0):
    """Shrink a figure-fraction rect by margin_mm on all four sides (3-mm whitespace)."""
    w_mm = fig.get_size_inches()[0] * 25.4
    h_mm = fig.get_size_inches()[1] * 25.4
    mx, my = margin_mm / w_mm, margin_mm / h_mm
    return [rect[0] + mx, rect[1] + my, rect[2] - 2 * mx, rect[3] - 2 * my]

def save_fig(fig, name, dpi=300):
    pdf, png, tif = FIGDIR/f"{name}.pdf", FIGDIR/f"{name}.png", FIGDIR/f"{name}.tif"
    fig.savefig(pdf, dpi=dpi, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=dpi)
    fig.savefig(tif, dpi=dpi, pil_kwargs={"compression": "tiff_lzw"})   # bitmap copy
    _png_to_rgb(png)
    _tif_to_rgb(tif)
    _strip_pdf_metadata(pdf)
    plt.close(fig)
    return pdf, png

def save_panel(fig, subdir, name, dpi=300):
    """Save a standalone panel (figures_panels/{subdir}/{name}.{ext}) — base name only."""
    base = FIGDIR.parent / "figures_panels" / subdir
    base.mkdir(parents=True, exist_ok=True)
    pdf, png, tif = base/f"{name}.pdf", base/f"{name}.png", base/f"{name}.tif"
    fig.savefig(pdf, dpi=dpi, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=dpi)
    fig.savefig(tif, dpi=dpi, pil_kwargs={"compression": "tiff_lzw"})
    _png_to_rgb(png)
    _tif_to_rgb(tif)
    _strip_pdf_metadata(pdf)
    plt.close(fig)
    return pdf, png

def _png_to_rgb(png):
    """Journal bitmap spec: RGB (matplotlib writes RGBA by default)."""
    from PIL import Image
    Image.open(png).convert("RGB").save(png, dpi=(300, 300))

def _tif_to_rgb(tif):
    """Journal bitmap spec: TIFF in RGB (matplotlib writes RGBA by default)."""
    from PIL import Image
    Image.open(tif).convert("RGB").save(tif, compression="tiff_lzw", dpi=(300, 300))

def _strip_pdf_metadata(pdf):
    """Remove Creator/Producer from the vector PDF (journal: no embedded metadata)."""
    try:
        import pikepdf
        with pikepdf.open(pdf, allow_overwriting_input=True) as p:
            with p.open_metadata(set_pikepdf_as_editor=False) as m:
                for k in list(m.keys()):
                    del m[k]
            p.save(pdf)
    except ImportError:
        pass   # pikepdf optional; matplotlib metadata already minimal

def panel_label(ax, letter, x=-0.24, dy_pt=18):
    """Panel letter left of the y-tick labels, 18 pt above the axes top (>= 6 mm
    clearance); 14 pt bold uppercase (user spec: >= si-hao/wu-hao size). Spanning
    panels at the figure's left edge get the letter at the axes' own top-left."""
    pos = ax.get_position()
    if pos.x0 + x * pos.width < 0.004:
        x = 0.0
    ax.annotate(letter, (x, 1.0), xycoords="axes fraction", xytext=(0, dy_pt),
                textcoords="offset points", fontsize=14, fontweight="bold",
                va="bottom", ha="left")

def pair_label(p):
    """'SRC->TGT' (or 'SRC -> TGT') with nephrology abbreviations, e.g. 'M2 macro→Monocyte'."""
    a, b = [x.strip() for x in p.replace(" -> ", "->").split("->")]
    return f"{CT_ABBR[a]}→{CT_ABBR[b]}"

def pval_text(p):
    return f"p = {p:.1e}" if p < 0.001 else f"p = {p:.4g}"

def sig_stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

def strip_pval_axis(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

def hline_legend(ax, handles, labels, loc="lower center", ncol=None, fontsize=MIN_FS, **kw):
    """Clean horizontal legend below a panel (dedicated space, never floating over data)."""
    if labels is None:
        labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, loc=loc, ncol=ncol or len(labels), fontsize=fontsize,
              handlelength=0.9, handletextpad=0.35, columnspacing=0.6,
              bbox_to_anchor=(0.5, -0.06), frameon=False, **kw)

def deoverlap_labels(ax, xs, ys, texts, fontsize=MIN_FS, pad_pt=1.5):
    """Greedy collision avoidance for direct labels; returns list of placed Text objects.
    Labels sorted by |x| (most extreme first) and nudged vertically to clear each other.
    Labels always sit to the RIGHT of their point so they never cross the y-axis labels."""
    import numpy as _np
    renderer = ax.figure.canvas.get_renderer()
    placed = []
    for x, y, t in sorted(zip(xs, ys, texts), key=lambda v: -abs(v[0])):
        txt = ax.annotate(t, (x, y), fontsize=fontsize, xytext=(4, 0), textcoords="offset points",
                          annotation_clip=False)
        bbox = txt.get_window_extent(renderer=renderer)
        for _ in range(60):
            hits = [p for p in placed if bbox.overlaps(p.get_window_extent(renderer=renderer))]
            if not hits:
                break
            txt.set_position((txt.get_position()[0], txt.get_position()[1] + bbox.height * 0.9))
            bbox = txt.get_window_extent(renderer=renderer)
        placed.append(txt)
    return placed
