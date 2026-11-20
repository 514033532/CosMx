"""Shared data loading for Fig 7 panels."""
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from common import SD

MANI = pd.read_csv(SD / "mif_images.csv", comment="#")
GDENS = pd.read_csv(SD / "mif_glom_densities.csv", comment="#")
FC = pd.read_csv(SD / "mif_foldchange.csv", comment="#")
MWU = pd.read_csv(SD / "mif_roi_mwu.csv", comment="#")
