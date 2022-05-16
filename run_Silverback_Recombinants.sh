#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

snakemake \
      -s Snakefile_Recombinants \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=240:00:00 -e logs -o logs" \
      --jobs 10 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 
