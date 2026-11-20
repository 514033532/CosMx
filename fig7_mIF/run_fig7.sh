#!/bin/bash
# one-script entry: fig7_mIF
cd "$(dirname "$0")"
PYTHONPATH="../lib:$PYTHONPATH" python3 plot_fig7.py "$@"
