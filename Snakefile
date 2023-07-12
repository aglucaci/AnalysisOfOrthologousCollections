"""``Snakemake`` file that runs the AOC alignment step.
Written by Alexander G Lucaci, 2023
"""
#----------------------------------------------------------------------------
# Imports
#----------------------------------------------------------------------------

import itertools
import os
import sys
import csv
import json
from pathlib import Path
from snakemake.utils import min_version
min_version("7.2")
from Bio import SeqIO
from ete3 import Tree
from Bio import Entrez
from ete3 import NCBITaxa
import pandas as pd
pd.options.mode.chained_assignment = None

#----------------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------------
configfile: 'config.yml'

print("# Loaded config yaml file")

with open("cluster.json", "r") as fh:
  cluster = json.load(fh)
  fh.close()
#end with

BASEDIR = os.getcwd()

print("# Loaded cluster json file")

Email = config["Email"]

Label = config["Label"]

Nucleotide_file = os.path.join(BASEDIR,
                               "data",
                               Label,
                               config["Nucleotide"])
                               
Protein_file = os.path.join(BASEDIR,
                            "data",
                            Label,
                            config["Protein"])
                            
CSV = os.path.join(BASEDIR,
                   "data",
                   Label,
                   config["CSV"])

print("# Using nucleotide data from:", Nucleotide_file)
print("# Using protein data from:", Protein_file)
print("# Using the analysis label:", Label)
print("# Using the metadata from:", CSV)

# Batch files
PREMSA  = os.path.join(BASEDIR,
                       "hyphy-analyses",
                       "codon-msa",
                       "pre-msa.bf")
                      
POSTMSA = os.path.join(BASEDIR,
                       "hyphy-analyses",
                       "codon-msa",
                       "post-msa.bf")

# Set output directory
print("# We are operating out of base directory:", BASEDIR)
       
OUTDIR_RESULTS = os.path.join(BASEDIR,
                              "results")
       
OUTDIR = os.path.join(BASEDIR,
                      "results",
                      Label)

print("# We will create and store results in:", OUTDIR)

# Create output dir.
os.makedirs(OUTDIR_RESULTS,
            exist_ok = True)
            
print("# Directory '% s' created" % OUTDIR_RESULTS)

os.makedirs(OUTDIR,
            exist_ok = True)
            
print("# Directory '% s' created" % OUTDIR)

# Set up configuration
PPN = cluster["__default__"]["ppn"]

# Hard-coded HyPhy settings
HYPHY = "hyphy"
HYPHYMPI = "HYPHYMPI"

CODON_OUTPUT = os.path.join(OUTDIR, Label)

CODONS_PY = os.path.join("scripts", "codons.py")

GARD_PARASE_NB = os.path.join("notebooks", "GARD_parse.ipynb")

#----------------------------------------------------------------------------
# Rule all
#----------------------------------------------------------------------------
def match_transcript_to_tree(TREE_NEWICK, accession):
    t = Tree(TREE_NEWICK, format=1)
    for leafname in t.get_leaf_names():
        if accession in leafname:
            return leafname
        #end if
    #end for
#end match

def ProcessLineages(transcript_accessions, DATA_DICT, TREE_NEWICK):
    count = 1
    for ACCESSION in transcript_accessions:
        skip = False
        
        for i in DATA_DICT.keys():
            if ACCESSION == DATA_DICT[i]["ACCESSION"]:
                skip = True
                break
            #end if
        #end for
        
        if skip == True:
            count += 1
            continue
        #end if
        
        try:
            handle = Entrez.esummary(db="nucleotide", id=ACCESSION, rettype="gb", retmode="text", retmax=1)
            records = Entrez.parse(handle)
        except Exception as e:
            print("# Error, sleeping", e)
            time.sleep(5)
            handle = Entrez.esummary(db="nucleotide", id=ACCESSION, rettype="gb", retmode="text", retmax=1)
            records = Entrez.parse(handle)
        #end try
        
        try:
            for record in records:
                TAXON_ID = record["TaxId"]
                print("#", count, "Processing transcript accession:", str(ACCESSION), "with NCBI Taxon ID:", str(TAXON_ID))
                ncbi = NCBITaxa()
                lineage = ncbi.get_lineage(TAXON_ID)
                names = ncbi.get_taxid_translator(lineage)
                leafname = ""
                leafname = match_transcript_to_tree(TREE_NEWICK, ACCESSION.replace(".", "_"))
                DATA_DICT[str(count)] = {"ACCESSION": ACCESSION, "TAXON_ID": TAXON_ID,
                                         "LINEAGE": [names[taxid] for taxid in lineage],
                                         "TITLE": record["Title"], "LEAFNAME": leafname}
                count += 1
            #end inner for
            handle.close
        except Exception as e:
            print("# Error (main):", ACCESSION, e, "\n")
        #end try
    #end outer for
    return DATA_DICT
#end method

def get_LineageColumn(lineages, loc):
    result = []
    for item in lineages:
        result.append(item[loc])
    #end for
    return result
#end method

#----------------------------------------------------------------------------
# Rule all
#----------------------------------------------------------------------------
rule all:
    input:
        CODON_OUTPUT,
        os.path.join(OUTDIR, Label + "_protein.fas"),
        os.path.join(OUTDIR, Label + "_nuc.fas"),
        os.path.join(OUTDIR, Label + "_protein.aln"),
        os.path.join(OUTDIR, Label + "_codons.fasta"),
        os.path.join(OUTDIR, Label + "_codons_duplicates.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta"),
        #os.path.join(OUTDIR, Label + "_codons.SA.fasta.treefile"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.GARD.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.best-gard"),
        os.path.join(OUTDIR, Label + ".1.codon.fas"),
        os.path.join(OUTDIR, Label + ".1.codon.fas.tree.nwk"),
        os.path.join(OUTDIR, Label + "_Annotated.csv")
#end rule all

print("# Moving on to processing rules")

#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------

rule get_codons:
    output:
        codons = CODON_OUTPUT
    params:
        Nuc = Nucleotide_file,
        Prot = Protein_file,
        Out = CODON_OUTPUT
    script:
        "scripts/codons.py"
#end rule

rule pre_msa:
    input: 
        codons = rules.get_codons.output.codons
    output: 
        protein_fas = os.path.join(OUTDIR, Label + "_protein.fas"),
        nucleotide_fas = os.path.join(OUTDIR, Label + "_nuc.fas")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} {PREMSA} --input {input.codons}"
#end rule 

rule mafft:
    input:
        protein = rules.pre_msa.output.protein_fas
    output:
        protein_aln = os.path.join(OUTDIR, Label + "_protein.aln")
    shell:
        "mafft --auto {input.protein} > {output.protein_aln}"
#end rule

rule post_msa:
    input: 
        protein_aln = rules.mafft.output.protein_aln,
        nucleotide_seqs = rules.pre_msa.output.nucleotide_fas  
    output: 
        codons_fas = os.path.join(OUTDIR, Label + "_codons.fasta"),
        duplicates_json = os.path.join(OUTDIR, Label + "_codons_duplicates.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} {POSTMSA} --protein-msa {input.protein_aln} --nucleotide-sequences {input.nucleotide_seqs} --output {output.codons_fas} --duplicates {output.duplicates_json}"
#end rule 

rule strike_ambigs:
   input:
       in_msa = rules.post_msa.output.codons_fas
   output:
       out_strike_ambigs = os.path.join(OUTDIR,
                                        Label + "_codons.SA.fasta")
   shell:
       "{HYPHY} scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule

#----------------------------------------------------------------------------
# IQ-TREE for ML tree inference
#----------------------------------------------------------------------------

rule iqtree: # Unrooted
    input:
        codons_fas = rules.strike_ambigs.output.out_strike_ambigs
    output:
        tree = os.path.join(OUTDIR, Label + "_codons.SA.fasta.treefile")
    shell:
        "iqtree -s {input.codons_fas} -T AUTO"
#end rule iqtree

#----------------------------------------------------------------------------
# Recombination detection
#----------------------------------------------------------------------------

rule recombination:
    input: 
        input = rules.strike_ambigs.output.out_strike_ambigs
    output: 
        output   = os.path.join(OUTDIR, Label + "_codons.SA.fasta.GARD.json"),
        bestgard = os.path.join(OUTDIR, Label + "_codons.SA.fasta.best-gard")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI}  GARD --alignment {input.input} --rv GDD --output {output.output}"
#end rule

#----------------------------------------------------------------------------
# Split out GARD partitions
#----------------------------------------------------------------------------

rule ParseNexus:
    input:
        input  = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR,
                              Label + ".1.codon.fas"),
        newick = os.path.join(OUTDIR,
                              Label + ".1.codon.fas.tree.nwk")
    params:
        output = os.path.join(OUTDIR, Label)
    run:
        # Load Partition coordinates
        """
        CHARSET span_1 = 1-2415;
        CHARSET span_2 = 2416-3003;
        """
        
        data = [line.strip() for line in open(input.input, "r").readlines() if "CHARSET" in line]
        coords = []
        for pos, i in enumerate(data):
            start = int(i.split(" ")[-1].split(";")[0].split("-")[0]) - 1
            stop  = int(i.split(" ")[-1].split(";")[0].split("-")[1]) - 1
            coords.append((start,stop))
        #end for
        print("# Number of partitions:", len(coords))
        print("# These are our coordinates:", coords)
        
        # Load Trees
        data  = [line.strip() for line in open(input.input, "r").readlines() if "TREE tree" in line]
        trees = []
        
        for n, tree in enumerate(data):
            # TREE tree_1 =
            newick = tree.split(" = ")[1]
            trees.append(newick)
            # Save tree to file
            t = Tree(newick, format=1)
            output_tree = params.output + ".".join(["." + str(n+1), "codon", "fas", "tree", "nwk"])
            t.write(format=1, outfile=output_tree)
        #end for
        
        # Report to user
        print("# Number of tree files (newick format):", len(trees))
        #print(trees)
        
        # Load records
        #records = SeqIO.parse(input.input, "nexus")
        
        #output_partitions
        for n, partition in enumerate(coords):
            print("# Working on partition", str(n+1), "with coordinates:", partition)
            # Load records
            records = SeqIO.parse(input.input, "nexus")
            output_records = []
            for record in records:
                output_record = record # Make a copy
                #print(len(output_record.seq))
                output_record.seq = output_record.seq[partition[0]: partition[1] + 1]
                output_records.append(output_record)
            #end for
            # Output records
            output_fasta = params.output + ".".join(["." + str(n + 1),
                                                     "codon",
                                                     "fas"])
            print("# Saving to FASTA file:", output_fasta)
            SeqIO.write(output_records, output_fasta, "fasta")
        #end for
                
#end rule

#----------------------------------------------------------------------------
# Selection analyses
#----------------------------------------------------------------------------

rule FEL:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.FEL.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} FEL --alignment {input.input} --output {output.output} --ci Yes"
#end rule

rule MEME:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.MEME.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} MEME --alignment {input.input} --output {output.output}"
#end rule

rule FUBAR:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.FUBAR.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} FUBAR --alignment {input.input} --output {output.output} "
#end rule

rule SLAC:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.SLAC.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} SLAC --alignment {input.input} --output {output.output}"
#end rule

rule FMM:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.FMM.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} FMM --alignment {input.input} --output {output.output}"
#end rule

rule BGM:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BGM.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BGM --alignment {input.input} --output {output.output}"
#end rule

rule aBSREL:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.aBSREL.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} ABSREL --alignment {input.input} --output {output.output}"
#end rule

rule aBSRELS:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.aBSRELS.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} ABSREL --alignment {input.input} --output {output.output} --srv Yes"
#end rule

rule BUSTED:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTED.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.input} --output {output.output} --srv No --starting-points 10"
#end rule

rule BUSTEDS:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTEDS.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.input} --output {output.output} --srv Yes --starting-points 10"
#end rule

rule BUSTEDMH:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTEDMH.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.input} --output {output.output} --srv No --starting-points 10 --multiple-hits Double+Triple"
#end rule

rule BUSTEDSMH:
    input:
        input = rules.recombination.output.bestgard
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTEDSMH.json")
    shell:
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.input} --output {output.output} --srv Yes --starting-points 10 --multiple-hits Double+Triple"
#end rule

#----------------------------------------------------------------------------
# Lineages
#----------------------------------------------------------------------------

rule GatherLineages:
    input:
        input_csv = CSV,
        tree      = rules.iqtree.output.tree
    params:
        email = Email,
        OUTDIR = OUTDIR
    output:
        output = os.path.join(OUTDIR, Label + "_Annotated.csv")
    run:
        Entrez.email = params.email
        df = pd.read_csv(input.input_csv)
        df.index += 1
        DATA_DICT = {}
        transcript_accessions = df['RefSeq Transcript accessions'].tolist()
        with open(input.tree, "r") as fh:
            TREE_NEWICK = fh.read()
        #end with
        DATA_DICT = ProcessLineages(transcript_accessions, DATA_DICT, TREE_NEWICK)
        df2 = pd.DataFrame.from_dict(DATA_DICT, orient="index")
        lineages = df2['LINEAGE'].tolist()
        #print(lineages)
        num_taxa = 5 # User-Set
        print("# Optimizing clade labels")

        count = 0
        BreakOut = False

        while True:
            lineageColumn = get_LineageColumn(lineages, count)
            lineageSet = set(lineageColumn)
            CountsLineageSet = []
            for animal in lineageSet:
                setCount = lineageColumn.count(animal)
                CountsLineageSet.append(setCount)
                #print(animal, setCount)
            #end for
            
            if min(CountsLineageSet) <= 5:
                break
            #end if
            
            count += 1
        #end while

        df2["CladeLabel"] = ""
        for index, row in df2.iterrows():
            #print(row['c1'], row['c2'])
            lineage = row['LINEAGE']
            for item in lineage:
                if item in lineageSet:
                    df2["CladeLabel"][index] = item
                #end if
            #end for
        #end for
                  
        print("# Saving Lineages:", output.output)
        df2.to_csv(output.output)
        
        # Output .clade files
        for index, row in df2.iterrows():
            leaf  = row["LEAFNAME"]
            CladeLabel = row["CladeLabel"]
            
            with open(os.path.join(OUTDIR, CladeLabel + ".clade"), "a") as fh:
                fh.write(str(leaf) + "\n")
                fh.close()
            #end with
        #end for
        
#end rule

#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------
