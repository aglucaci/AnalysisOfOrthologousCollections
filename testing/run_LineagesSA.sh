#!/bin/bash
#PBS -l walltime=999:00:00

#@Author: Alexander G Lucaci
#@Usage qsub -V -l nodes=1:ppn=16 -q epyc run_LineagesSA.sh

WD=/home/aglucaci/AOC/results/BDNF/Recombinants
ALN=/home/aglucaci/AOC/results/BDNF/Recombinants/BDNF_codons_RDP_recombinationFree.fas

HYPHY=/home/aglucaci/hyphy-develop/HYPHYMPI
RES=/home/aglucaci/hyphy-develop/res

for tree in $WD/*.treefile.*; do 
    tbase="$(basename -- $tree)"
    echo "# basename: "$tbase
    #label=${tbase%BDNF_codons_RDP_recombinationFree.fas.treefile}
    label=${tbase/BDNF_codons_RDP_recombinationFree.fas.treefile./}

    echo "# Label: "$label

    # Do MEME with branches labelled
    # -----------------------------------------------------------------
    outputMEME="$ALN"."$label".MEME.json
    echo "# Output file: "$outputMEME

    if [ ! -f $outputMEME ]; then
        echo "# Execute hyphy command"
        echo mpirun -np 16 $HYPHY LIBPATH=$RES MEME --alignment $ALN --tree $tree --branches $label --output $outputMEME
        mpirun -np 16 $HYPHY LIBPATH=$RES MEME --alignment $ALN --tree $tree --branches $label --output $outputMEME
    else
        echo "# output file exists, verify results"
    fi
  
    # Debug
    #continue

    # Do FEL with branches labelled
    # -----------------------------------------------------------------
    outputFEL="$ALN"."$label".FEL.json
    echo "# Output file: "$outputFEL

    if [ ! -f $outputFEL ]; then
        echo "# Execute hyphy command"
        echo mpirun -np 16 $HYPHY LIBPATH=$RES FEL --alignment $ALN --tree $tree --branches $label --output $outputFEL
        mpirun -np 16 $HYPHY LIBPATH=$RES FEL --alignment $ALN --tree $tree --branches $label --output $outputFEL

    else
        echo "# output file exists, verify results"
    fi

    # Do RELAX with branches labelled
    # -----------------------------------------------------------------
    outputRELAX="$ALN"."$label".RELAX.json
    echo "# Output file: "$outputRELAX
    #echo ""
    if [ ! -f $outputRELAX ]; then
        echo "# Execute hyphy command"

        echo mpirun -np 16 $HYPHY LIBPATH=$RES RELAX --alignment $ALN --tree $tree --test $label --mode "Group mode" --output $outputRELAX
        mpirun -np 16 $HYPHY LIBPATH=$RES RELAX --alignment $ALN --tree $tree --test $label --mode "Group mode" --output $outputRELAX
    else
        echo "# output file exists, verify results"
    fi
  
    # Do CFEL with branches labelled
    # -----------------------------------------------------------------
    outputCFEL="$ALN"."$label".CFEL.json
    echo "# Output file: "$outputCFEL
    #echo ""
    if [ ! -f $outputCFEL ]; then
        echo "# Execute hyphy command"
        echo mpirun -np 16 $HYPHY LIBPATH=$RES CONTRAST-FEL --alignment $ALN --tree $tree --branch-set $label --output $outputCFEL
        mpirun -np 16 $HYPHY LIBPATH=$RES CONTRAST-FEL --alignment $ALN --tree $tree --branch-set $label --output $outputCFEL
    else
        echo "# output file exists, verify results"
    fi
    echo "" 
done





exit 0
