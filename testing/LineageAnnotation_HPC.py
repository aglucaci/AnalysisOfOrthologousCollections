#!/usr/bin/env python
# coding: utf-8

# In[1]:

# Imports
import sys, os
from Bio import Entrez
import json
from ete3 import NCBITaxa
import pandas as pd
from ete3 import Tree


# In[2]:


# Declares
CSV_FILE = "../data/BDNF/BDNF_orthologs.csv"

DATA_DICT = {}

# The following looks up NCBI Taxon ID from Accession IDs
# Then uses that Taxon ID to get the complete taxnomic lineage of of that species.
Entrez.email = "aglucaci@gmail.com"
#Entrez.email = snakemake.params.Entrez_Email

# In[4]:
TREE_FILE = "../results/BDNF/Recombinants/BDNF_codons_RDP_recombinationFree.fas.treefile"

# Read file, get newick
with open(TREE_FILE, "r") as fh:
    TREE_NEWICK = fh.read()
#end with

# In[5]:

df = pd.read_csv(CSV_FILE)
df.index += 1

# In[6]:

transcript_accessions = df['RefSeq Transcript accessions'].tolist()

#len(transcript_accessions)

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
    


# In[8]:


def main(transcript_accessions):
    #DATA_DICT = {}
    global DATA_DICT, TREE_NEWICK
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
            print("# Error (in main):", ACCESSION, e, "\n")
            #print(names)
            #print(e)
            #main(TREE_NEWICK) # Recursive.
            #break
        #end try
        
    #end outer for
    
    #return DATA_DICT
#end method


# In[9]:


main(transcript_accessions)


# In[10]:


# Check for errors here
print(len(transcript_accessions)) # Total
print(len(DATA_DICT.keys())) # What we have identified
# They need to match


# In[11]:


df2 = pd.DataFrame.from_dict(DATA_DICT, orient="index")
#df2


# In[32]:


lineages = df2['LINEAGE'].tolist()
num_taxa = 20

for i in range(len(lineages[0])):
    
    to_add = []
    
    for species in lineages:
        to_add.append(species[i])
    #end inner for
    
    df2["Annotation"] = to_add
    
    a = df2["Annotation"].nunique() # number of unique taxa.
    
    if a >= num_taxa: break
        
    print("#",i, "=", df2["Annotation"].nunique())
#end outer for



for item in set(df2["Annotation"].to_list()):
    #print(item)
    df3 = df2[df2["Annotation"] == item]
    #print(df3)
    df3 = df3.dropna()
    
    #if len(df3["ACCESSION"].to_list()) < 3: 
    #    continue
    #print(df3["LEAFNAME"].to_list())
    #df3.to_csv(item + ".csv", index=True)
    
    output_file_clade = os.path.join(*["..", "results", "BDNF", "Recombinants", item.replace(" ", "_") + ".clade"])
    
    #print("# Saving", str(len(df3["LEAFNAME"].to_list())), "to:", item + ".clade")
    print("# Saving", str(len(df3["LEAFNAME"].to_list())), "to:", output_file_clade)
    
    with open(output_file_clade, "w") as fh:
        fh.write("\n".join(df3["LEAFNAME"].to_list()))
        #for row in df3["LEAFNAME"].to_list():
        #    fh.write(row)
    #end with
#end for


# END OF FILE
