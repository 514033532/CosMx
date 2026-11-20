#!/bin/bash
# one-script entry: fig6_communication
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_fig6.py "$@"
