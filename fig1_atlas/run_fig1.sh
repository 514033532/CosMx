#!/bin/bash
# one-script entry: fig1_atlas
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_fig1.py "$@"
