#!/usr/bin/env zsh --login

# Script to activate python environment (for data analysis) and
# run perl script (to build and call C++ project).
#
# Command line arguments:
#     $0: ./run.sh
#     $1:  filename of perl script to call (no extension or prefix)
#
# EXAMPLE: ./run.sh collinear-swimmer-wall

# exit immediately if a command exits with a non-zero status.
set -e

# activate conda environment
conda activate bodies-in-potential-flow

# run script in shell
exec perl scripts/"$@".pl