#!/bin/bash

clear
echo "alpha v0.1 --- "
echo "2022 - Analysis of Orthologous Collections (AOC)."
echo "@Author: Alexander G. Lucaci"
echo ""

# Set up the pipeline failure expectations.
set -euo pipefail

# Downloading hyphy-analyses.
FOLDER="hyphy-analyses"
URL="https://github.com/veg/hyphy-analyses.git"

if [ ! -d "$FOLDER" ] ; then
    git clone "$URL" "$FOLDER"
fi

echo "Initialized --- "

# Uncomment if you want to generate an analysis DAG file.
#snakemake --forceall --dag | dot -Tpdf > AOC_Alignment.pdf

echo "Creating 'log' directory"
mkdir -p logs

echo "Executing LOCAL Snakemake command"

# Snakemake command ----
snakemake \
      -s Snakefile_Alignment \
      --jobs 2 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 

# End of file  -- done.
