#!/bin/bash
# one-script entry: fig5_3D_DE
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_fig5.py "$@"
