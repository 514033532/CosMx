"""Snapshot the current figure set into a dated round folder.

Usage: python3 snapshot_round.py "<one-line summary of this correction round>"
Creates figures/rounds/YYYYMMDD_HHMM/ containing all current figures (pdf/png/tif),
standalone panels, and a README with the round summary. Run once per correction round.
"""
import sys
import shutil
from pathlib import Path
from datetime import datetime
import common

def main():
    summary = sys.argv[1] if len(sys.argv) > 1 else "(no summary)"
    stamp = datetime.now().strftime("%Y%m%d_%H%M")
    round_dir = common.FIGDIR / "rounds" / stamp
    round_dir.mkdir(parents=True, exist_ok=True)
    for f in sorted(common.FIGDIR.glob("*")):
        if f.is_file() and f.suffix.lower() in (".pdf", ".png", ".tif"):
            shutil.copy(f, round_dir / f.name)
    panels_root = common.FIGDIR.parent / "figures_panels"
    if panels_root.exists():
        for d in sorted(panels_root.glob("fig*")):
            if d.is_dir():
                dst = round_dir / f"panels_{d.name}"
                dst.mkdir(exist_ok=True)
                for f in sorted(d.glob("*")):
                    if f.is_file() and f.suffix.lower() in (".pdf", ".png", ".tif"):
                        shutil.copy(f, dst / f.name)
    (round_dir / "README.md").write_text(
        f"# 轮次归档 {stamp}\n\n本轮回纠内容：{summary}\n\n"
        f"包含：全部当前图件（fig1-5 人工版、fig6/7/SFig1/2 重绘版、fig8 草稿）+ 独立面板。\n")
    print(f"round snapshot -> {round_dir}")

if __name__ == "__main__":
    main()
