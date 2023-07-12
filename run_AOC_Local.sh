#!/bin/bash

set -euo pipefail

# Downloading hyphy-analyses.
FOLDER="hyphy-analyses"

URL="https://github.com/veg/hyphy-analyses.git"

if [ ! -d "$FOLDER" ] ; then
    git clone "$URL" "$FOLDER"
fi

# Set up logs directory
LOGS="logs"

if [ ! -d "$LOGS" ] ; then
    mkdir -p $LOGS
fi

printf "# Running AOC Snakemake pipeline \n"

# Initial phase, quality control, alignment, recombination-detection.
snakemake \
      -s Snakefile \
      --cluster-config cluster.json \
      --jobs 1 all \
      --rerun-incomplete \
      --cores all \
      --keep-going \
      --reason \
      --latency-wait 120
      
#exit 0

# Run Selection Analyses on recombination-free files
snakemake \
      -s Snakefile_Recombinants \
      --cluster-config cluster.json \
      --jobs 1 all \
      --rerun-incomplete \
      --cores all \
      --keep-going \
      --reason \
      --latency-wait 120

exit 0
# End of file

