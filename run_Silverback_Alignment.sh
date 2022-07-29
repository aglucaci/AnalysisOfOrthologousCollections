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

echo "Executing HPC Snakemake command"

snakemake \
      -s Snakefile_Alignment \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=720:00:00 -e logs -o logs" \
      --jobs 10 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 

# End of file 
