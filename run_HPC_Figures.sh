#!/bin/bash

echo "Analysis of Orthologous Collections (AOC)"
echo "@Author: Alexander G. Lucaci"
echo ""

# Set up the pipeline failure expectations.
set -euo pipefail

echo "Executing HPC Snakemake command, Generating figures and tables"

snakemake \
      -s Snakefile_Figures \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 8 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 120 

exit 0 

# End of file 
