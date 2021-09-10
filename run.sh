#!/usr/bin/env bash --login
set -e

# activate conda environment and let the following process take over
conda activate bodies-in-potential-flow
exec scripts/"$@".pl