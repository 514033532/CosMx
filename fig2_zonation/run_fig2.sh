#!/bin/bash
# one-script entry: fig2_zonation
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_fig2.py "$@"
