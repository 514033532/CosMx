"""Shared data loading for Fig 6 panels (one source of truth per figure)."""
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common
from common import SD

CS = pd.read_csv(SD / "comm_strength.csv", comment="#")
LN = pd.read_csv(SD / "lr_num_de.csv", comment="#")
LV = pd.read_csv(SD / "lr_value_de.csv", comment="#")
CD = pd.read_csv(SD / "centrality_de.csv", comment="#")
LS = pd.read_csv(SD / "lr_pair_samples.csv", comment="#")
CMP = "Group2_DN-VS-NC"
