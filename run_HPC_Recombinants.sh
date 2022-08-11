#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

snakemake \
      -s Snakefile_Recombinants \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 8 all \
      --rerun-incomplete \
      --keep-going \
      --use-conda \
      --reason \
      --latency-wait 60 
