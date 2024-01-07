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
# Run Selection Analyses on recombination-free files
snakemake \
      -s Snakefile \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 4 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 300
      
# Run Selection Analyses on recombination-free files
snakemake \
      -s Snakefile_Recombinants \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 4 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 300

exit 0
# End of file

