#!/bin/bash
clear

#set -euo pipefail

# Downloading hyphy-analyses.
FOLDER="software/hyphy-analyses"
URL="https://github.com/veg/hyphy-analyses.git"

if [ ! -d "$FOLDER" ] ; then
    git clone "$URL" "$FOLDER"
fi

# Set up logs directory
LOGS="logs"

if [ ! -d "$LOGS" ] ; then
    mkdir -p $LOGS
fi

#VERSION=0.1.0

###############################################################################
printf "###############################################################################\n"
printf "# Running the AOC Snakemake pipeline \n"
printf "###############################################################################\n"
###############################################################################

# Initial phase, quality control, alignment, recombination-detection.
# Run Selection Analyses on recombination-free files
snakemake \
      -s workflow/Snakefile \
      --cluster-config config/cluster.json \
      --jobs 2 all \
      --cores all \
      --keep-going \
      --reason \
      --latency-wait 60 \
      --rerun-incomplete

echo ""

###############################################################################
printf "###############################################################################\n"
printf "# Running the AOC Snakemake pipeline - Recombination-free analysis \n"
printf "###############################################################################\n"
###############################################################################
      
# Run Selection Analyses on recombination-free files
snakemake \
      -s workflow/Snakefile_SelectionAnalysis \
      --cluster-config config/cluster.json \
      --jobs 4 all \
      --cores all \
      --keep-going \
      --reason \
      --latency-wait 60 \
      --rerun-incomplete

echo ""

###############################################################################
printf "###############################################################################\n"
printf "# Running the AOC Snakemake pipeline - Visualization and Summary \n"
printf "###############################################################################\n"
###############################################################################

# Generate all output tables and charts.
snakemake \
      -s workflow/Snakefile_SummarizeResults \
      --cluster-config config/cluster.json \
      --jobs 1 all \
      --cores all \
      --keep-going \
      --reason \
      --latency-wait 60

exit 0

###############################################################################
# End of file
###############################################################################
