#!/bin/bash


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

###############################################################################
printf "###############################################################################\n"
printf "# Running the AOC Snakemake pipeline \n"
printf "###############################################################################\n"
###############################################################################

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
    
###############################################################################
printf "###############################################################################\n"
printf "# Running the AOC Snakemake pipeline - Recombination-free analysis \n"
printf "###############################################################################\n"
###############################################################################
  
# Run Selection Analyses on recombination-free files
snakemake \
      -s workflow/Snakefile_SelectionAnalysis \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 4 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 300

###############################################################################
printf "###############################################################################\n"
printf "# Running the AOC Snakemake pipeline - Visualization and Summary \n"
printf "###############################################################################\n"
###############################################################################

# Generate all output tables and charts.
snakemake \
      -s workflow/Snakefile_SummarizeResults \
      --cluster-config config/cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 1 all \
      --cores all \
      --keep-going \
      --reason \
      --latency-wait 60

exit 0

###############################################################################
# End of file
###############################################################################
