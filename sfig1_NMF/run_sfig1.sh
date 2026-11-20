#!/bin/bash
# one-script entry: sfig1_NMF
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_sfig1.py "$@"
