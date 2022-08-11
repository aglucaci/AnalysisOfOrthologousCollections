#!/usr/bin/env python
# coding: utf-8

# Imports
import sys, os
from Bio import Entrez
import json
from ete3 import NCBITaxa
import pandas as pd
from ete3 import Tree
#from tqdm import tqdm
#import urllib2
import time

# Declares
RESULTS_DIR     = sys.argv[1]
CSV_FILE        = sys.argv[2]
TREE_FILE       = sys.argv[3]
TREE_FILE_LOG   = TREE_FILE + ".log"

# The following looks up NCBI Taxon ID from Accession IDs
# Then uses that Taxon ID to get the complete taxnomic lineage of of that species.
Entrez.email = "alexander.lucaci@temple.edu"

# Read file, get newick
if os.stat(TREE_FILE).st_size == 0:
    print("# Empty tree file:", TREE_FILE)
    sys.exit(1)
#end if

with open(TREE_FILE, "r") as fh:
    TREE_NEWICK = fh.read()
#end with

print("# Tree loaded:", TREE_FILE)
DATA_DICT = {}

if os.stat(CSV_FILE).st_size == 0:
    print("# Empty orthologs (csv) file:", TREE_FILE)
    sys.exit(1)

print("# Reading orthologs (csv) file:", CSV_FILE)
df = pd.read_csv(CSV_FILE)
df.index += 1
# Process
transcript_accessions = df['RefSeq Transcript accessions'].tolist()
# In[7]:
# Helper function
def match_transcript_to_tree(TREE_NEWICK, accession):
    t = Tree(TREE_NEWICK, format=1)
    # uncomment to debug
    #print("# in 'match_transcript_to_tree'", accession)
    for leafname in t.get_leaf_names():
        if accession in leafname:
            # uncomment to debug
            #print("# Returning:", leafname)
            return leafname
        #end if
    #end for
#end match


def main(transcript_accessions, DATA_DICT, TREE_NEWICK):
    count = 1
    for ACCESSION in transcript_accessions:
        #print("# Checking ACCESION:", ACCESSION)
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
                # uncomment to debug
                #print(lineage) # returns a list of taxon id's
                names = ncbi.get_taxid_translator(lineage)
                # uncomment to debug
                #print(names) # translates the taxon id's into clade names
                # uncomment to debug
                #print("# ACCESSION:", ACCESSION)
                leafname = ""
                leafname = match_transcript_to_tree(TREE_NEWICK, ACCESSION.replace(".", "_"))
                # uncomment to debug
                #print(leafname)
                DATA_DICT[str(count)] = {"ACCESSION": ACCESSION, "TAXON_ID": TAXON_ID, 
                                         "LINEAGE": [names[taxid] for taxid in lineage], 
                                         "TITLE":record["Title"], "LEAFNAME": leafname}
                count += 1
            #end inner for
            handle.close
        except Exception as e:
            print("# Error (main):", ACCESSION, e, "\n")
            #print(names)
            #print(e)
            #main(TREE_NEWICK) # Recursive.
            #break
        #end try
    #end outer for
    #return DATA_DICT
#end method
# In[11]:
# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

main(transcript_accessions, DATA_DICT, TREE_NEWICK)

# ---------------------------------------------------------------------
# ## Second phase
# ---------------------------------------------------------------------
df2 = pd.DataFrame.from_dict(DATA_DICT, orient="index")
lineages = df2['LINEAGE'].tolist()
num_taxa = 20 # User-Set

for i in range(len(lineages[0])):
    to_add = []
    for species in lineages:
        to_add.append(species[i])
    #end inner for
    df2["Annotation"] = to_add
    a = df2["Annotation"].nunique()
    if a >= num_taxa: break
    print("#",i, "=", df2["Annotation"].nunique())
#end outer for
# In[15]:

def optimize_lineages(df2, seed=5):
    # no lineage should be >= 40% of species.
    lineages = df2['LINEAGE'].tolist()
    num_taxa = seed
    print("## Seed =", seed)
    for i in range(len(lineages[0])):
        to_add = []
        for species in lineages:
            to_add.append(species[i])
        #end inner for
        df2["Annotation"] = to_add
        a = df2["Annotation"].nunique()
        if a >= num_taxa: 
            break
        #print("#",i, "=", df2["Annotation"].nunique())
    #end outer for
    # Do a check here.
    threshold = 0.4
    total = len(df2["Annotation"].tolist())
    for item in df2.groupby('Annotation').Annotation.count():
        calc = item / total
        print(calc)
        if calc > threshold:
            print(item / total, "exceeded threshold") 
            seed += 1
            optimize_lineages(df2, seed)
        else:
            pass
        #end if
    #end for
    # if it fails, make this recursive.
    # If we make it here, success
    print("# final seed =", seed)
    return df2
#end method

df2 = optimize_lineages(df2)
threshold = 0.4
total = len(df2["Annotation"].tolist())

for item in df2.groupby('Annotation').Annotation.count():
    calc = item / total
    print(item / total)
    if calc > threshold:
        print(item / total, "exceeded threshold") 


json_data = {}
for item in set(df2["Annotation"].to_list()):
    json_data[item] = {}

for item in set(df2["Annotation"].to_list()):
    df3 = df2[df2["Annotation"] == item]
    df3 = df3.dropna()
    #print("# Saving", str(len(df3["LEAFNAME"].to_list())), "to:", item + ".clade")
    # RESULTS_DIR
    out_clade_file = os.path.join(RESULTS_DIR, item.replace(" ", "_") + ".clade")
    if len(df3["LEAFNAME"].to_list()) <= 2:
        continue
    print("# Saving", str(len(df3["LEAFNAME"].to_list())), "to:", out_clade_file)
    with open(out_clade_file, "w") as fh:
        fh.write("\n".join(df3["LEAFNAME"].to_list()))
    with open(TREE_FILE_LOG, "a") as fh2:
        fh2.write("\n".join(df3["LEAFNAME"].to_list()))

    json_data[item] = df3["LEAFNAME"].to_list()
#end for

out_clade_file_json = os.path.join(RESULTS_DIR, "CladeLabels.json")
print("# Writing JSON file to:", out_clade_file_json)
with open(out_clade_file_json, 'w') as outfile:
    #outfile.write(json_data)
    json.dump(json_data, outfile)






