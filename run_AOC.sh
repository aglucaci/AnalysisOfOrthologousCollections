#!/bin/bash

# Try to set up the pipeline cleanly.

# Step 1 (Alignment + Recombination detection)
bash run_Silverback_Alignment.sh

# Step 2 (Selection analyses)
bash run_Silverback_Recombinants.sh

# Step 3 (Lineage annotation, and selection analyses)
bash run_Silverback_Recombinants_Lineages.sh

# Step 4 (Generate figures + figure legends)
bash run_Silverback_Figures.sh


exit 0
