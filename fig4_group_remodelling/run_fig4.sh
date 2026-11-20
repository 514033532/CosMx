#!/bin/bash
# one-script entry: fig4_group_remodelling
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_fig4.py "$@"
