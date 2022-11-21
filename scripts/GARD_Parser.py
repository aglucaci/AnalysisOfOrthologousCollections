#!/usr/bin/env python
# coding: utf-8

# ## Imports 

# In[1]:

from Bio import SeqIO
from Bio.Seq import Seq
import json
import os

# ## Declares

# In[2]:

# Empty

# ## Helper functions

# In[3]:


# Helper functions
def gard_parser(combo, best_gard, codon_MSA_in):    
    global BASEDIR
    data = [d.strip() for d in open(best_gard).readlines() if "CHARSET" in d]
    coords = []
    if len(data) > 1:
        for pos, i in enumerate(data):
            ## hyphy coordinates are 1 indexed and INCLUSIVE ##
            ## so to get to "true" python coords you need to -1 to start and -1 to stop ##
            ## then you need to +1 to python "stop" coord to get the actual length ## 
            #hyphy_start = i.split(' ')[-1].split(';')[0].split('-')[0]
            #hyphy_stop = i.split(" ")[-1].split(";")[0].split("-")[1]
            #print(f"HYPHY coords for partition {pos +1} --> start {hyphy_start} | stop {hyphy_stop}")
            ## in python (start-1) | (end-1) to get the true coords ##
            ## to get the correct python index[X:Y] (zero index, EXLUSIVE) need to do [start-1,stop] ##
            start = int(i.split(" ")[-1].split(";")[0].split("-")[0])-1
            stop = int(i.split(" ")[-1].split(";")[0].split("-")[1])-1
            coords.append((start,stop))        
        #end for
    ## print out fasta file
    else:
        print(data)
    #end if
    ## now use the coords to pull out seqs in the codon aware MSA ##
    ## write a json for partition and indices associated ## 
    index_data = {}
    for pos, c in enumerate(coords):
        temp_dict = {}
        ## want list of 3 ##
        n = 3
        recs = list(rec for rec in SeqIO.parse(codon_MSA_in, "fasta"))      
        ## actual codons ##
        #codon_map = [list(recs[0].seq)[i*n:(i+1) *n] for i in range((len(list(recs[0].seq)) + n - 1) // n )]
        index_map = [list(range(len(recs[0].seq)))[i*n:(i+1) *n] for i in range((len(list(recs[0].seq)) + n - 1) // n )]
        #print(index_map)
        old_start = c[0]
        old_stop = c[1]
        new_start = ''
        new_stop = ''
        ## getting new start ##
        for imap in index_map:
            if old_start in imap:
                for p, j in enumerate(imap):
                    if old_start == j:
                        if p !=0:
                            new_start = imap[2]+1
                        else:
                            new_start = imap[0]
                        #end if
                    #end if
                #end for
            elif old_stop in imap:
                new_stop = imap[0]
            else:
                continue 
            #end if
        ## sanity check 
        print(f"# Sanity check: OLD {old_start}, {old_stop} | NEW {new_start}, {new_stop} | NEW div by 3 {(int(new_stop)-int(new_start))/3}")
        #print(f"OLD {old_start}, {old_stop} | NEW {new_start}, {new_stop}")
        temp_dict[pos+1] = [i+1 for i in list(range(int(new_start),int(new_stop)))]
        #print(temp_dict)
        index_data.update(temp_dict)
        
        ### ~ CHANGE TO PLACE YOU WANT OUTPUT ~ ###
        codon_out = os.path.join(BASEDIR, "results", label, "%s.%s.codon.fas" % (combo, str(pos+1)))
        
        print("# Saving partition to:", codon_out)
        with open(codon_out, "w") as out_f:
            ## just grabbing one record right now ## 
            #for record in recs[0:1]:
            for record in recs:
                ## now grabbing actual partition ##
                partition = record[int(new_start):int(new_stop)] + Seq("---")
                ## double sanity check ##
                #print(len(partition),len(partition)/3)
                #print(partition.seq)

                if set(partition.seq) == set("-"): #100% gap sequence
                    continue
                #end if

                out_f.write(">{}\n{}\n".format(partition.id, partition.seq))
            #end for
        #end with
        
    #print(index_data)
    ## write json here ## 
    #json_title = "CATS_8-2021/CATS_8-20.codons/FCOV1_%s.PARTITION_info.json" % (combo)
    #with open(json_title, 'w') as json_file:
    #      json.dump(index_data, json_file, indent=4)
    #temp = {}
    #temp[combo] = list(index_data.keys())
    ##print(temp)
    print(f"\nDone with {combo}\n")
    #return temp


# # Main subroutine

# In[4]:

# Input from the snakemake pipeline.
BASEDIR = snakemake.params.BASEDIR
label =  snakemake.params.genelabel

# Grab the best-gard and input-msa, we need these.
# We are going from the best-gard charsets,
# And restoring the full-fasta file with breakpoints
# As were determined from the downsampled fasta.
# Test on mammalian_REV3L_codons.SA.fasta.cluster.fasta.best-gard
best_gard = os.path.join(BASEDIR, "results", label, label + "_codons.SA.fasta.cluster.fasta.best-gard")
input_msa = os.path.join(BASEDIR, "results", label, label + "_codons.SA.fasta") # Whole MSA not the clustered version

# Print to user
print("# Processing best-gard:", best_gard)
print("# Input Alignment:", input_msa)

# Do work.
gard_parser(label, best_gard, input_msa)

# ## End of file
