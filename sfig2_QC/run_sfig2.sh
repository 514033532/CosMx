#!/bin/bash
# one-script entry: sfig2_QC
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_sfig2.py "$@"
