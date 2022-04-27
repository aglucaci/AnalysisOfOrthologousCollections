# Imports

import os
import sys
import argparse
import json
import shutil
import random


arguments = argparse.ArgumentParser(description='Cluster an MSA with genetic distance (TN93)')

arguments.add_argument('-i', '--input',            help = 'MSA file to process',                                  required = True, type = str )
arguments.add_argument('-o', '--output_fasta',           help = 'Output json file',                                     required = True, type = str)

settings = arguments.parse_args()


def cluster_to_fasta (in_file, out_file, ref_seq = None):
    # assert that in_file exists, otherwise this will crash if tn93-cluster did not run.

    with open (in_file, "r") as fh:
        cluster_json = json.load (fh)
        #print (colored('Running ... converting representative clusters to .FASTA\n', 'cyan'))
        check_uniq = set ()
        print("# Saving to fasta:", out_file)
        with open (out_file, "w") as fh2:
            for c in cluster_json:
                cc = c['centroid'].split ('\n')
                if ref_seq:
                    if ref_seq in c['members']:
                        cc[0] = ">" + ref_seq
                        print ("\n".join (cc), file = fh2)
                        continue
                    #end if
                #end if
                seq_id = cc[0]
                while seq_id in check_uniq:
                        seq_id = cc[0] + '_' + ''.join(random.choices ('0123456789abcdef', k = 10))
                #end while
                check_uniq.add (seq_id)
                #print (seq_id,"\n",cc[1], file = fh2)
                print (seq_id + "\n" + cc[1].replace(" ", "") + "\n", file = fh2)
            #end for
         #end with
    #end with
    return (os.path.getmtime(out_file), len(cluster_json))
#end method


input_file = settings.input
compressed_fasta = settings.output_fasta 

cluster_json = input_file

#input_stamp, cluster_count = cluster_to_fasta (cluster_json, compressed_fasta)
input_stamp, cluster_count = cluster_to_fasta (cluster_json, compressed_fasta)


# Exit
