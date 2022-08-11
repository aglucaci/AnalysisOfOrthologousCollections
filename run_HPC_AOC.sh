#!/bin/bash

# Try to set up the pipeline cleanly.
# Move hyphy analyses step here

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


# Step 1 (Alignment + Recombination detection)
bash run_HPC_Alignment.sh

# Step 2 (Selection analyses)
bash run_HPC_Recombinants.sh

# Step 3 (Lineage annotation, and selection analyses)
bash run_HPC_Recombinants_Lineages.sh

# Step 4 (Generate figures + figure legends)
bash run_HPC_Figures.sh


exit 0





