#!/bin/bash
# one-script entry: fig3_zone_pathology
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_fig3.py "$@"
