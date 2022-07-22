# Analysis Of Orthologous Collections (AOC)

```alpha v.0.03 last updated 07.21.2022```

This repository expands on work oringally done for https://github.com/aglucaci/GeneInvestigator. 

A simple application to interrogate the molecular evolution of a single gene. AOC allows for the inclusion of recombination detection, a powerful force in shaping gene evolution and interpreting analytic results. As well, it allows for lineage assignment and annotation. This feature (lineage assignment) allows between group comparison of selective pressures.

## Installation and dependencies
This application is currently designed to run in an HPC environment due to the computational cost of selection analyse.

There is an assumption that the freely available Anaconda software is installed on your machine.

You will also need to download the standalone hyphy-analyses repository (https://github.com/veg/hyphy-analyses). 

Make sure to modify the config.yml file to point to the correct directory on your system

### To install -- Steps necessary to complete before running
1. `git clone https://github.com/aglucaci/AnalysisOfOrthologousCollections.git AOC`
2. `cd AOC`
3. `conda env create -f environment.yml`. This will create a virtual environment called (GeneInvestigator) with the necessary dependencies.
4. At this point, run `conda activate AOC` and your environment will be ready to go.

## Data retrieval via NCBI Orthologs
Here, we rely on the NCBI Ortholog database. 

For example, if we are interested in the TP53 gene: https://www.ncbi.nlm.nih.gov/gene/7157/ortholog/?scope=117570&term=TP53

Download all information: Tabular data, RefSeq Transcripts, and RefSeq Protein. 

This is typically done as one gene per species, but all transcripts per species is also available.

## Pipeline

1. Codon-aware alignment from protein and gene transcript files
2. Recombination detection (via RDP)
3. Tree inference (ML Tree inference, preferably IQ-TREE or RAxML).
4. Selection analyses: this includes (MEME, FEL, FUBAR, BUSTEDS MEME, aBSREL, SLAC, BGM, PRIME, FMM)
5. Lineage assignment and tree annotation (Custom script)
6. Selection analyses on lineages (RELAX, CFEL)
7. Visualization

## Results
The following are JSON files produced by HyPhy analyses. 

These can be visualized by the appropriate module from HyPhy Vision (http://vision.hyphy.org/). 

Analysis file names contain the method used (SLAC, FEL, PRIME, FADE, MEME, CFEL, etc), and if appropriate -- the set of branches to which the analysis was applied.

Using the TP53 example above.

```
── results/TP53/
│   ├── TP53.FEL.json
│   ├── TP53.FUBAR.json
│   ├── TP53.BUSTEDS.json
│   ├── TP53.MEME.json
│   ├── TP53.ABSREL.json
│   ├── TP53.SLAC.json
│   ├── TP53.BGM.json
│   ├── TP53.PRIME.json
│   ├── TP53.FMM.json
```

## Removing the AOC environment from conda
```conda env remove --name AOC```

## Starting a tmux session for AOC
```tmux new -s AOC ```
