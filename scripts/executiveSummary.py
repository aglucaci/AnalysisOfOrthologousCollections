"""
Summary script for HyPhy JSON files from various methods.

@Author: Alexander G. Lucaci

"""

# Imports
import pandas as pd
import numpy as np
import os
import json
import glob
from os.path import exists
import sys
import argparse
import statsmodels
import statsmodels.api

# Declares

jsonFileFEL   = snakemake.input.inputFEL
jsonFileFUBAR = snakemake.input.inputFUBAR
jsonFileMEME  = snakemake.input.inputMEME
#jsonFileBUSTEDS  = snakemake.input.inputBUSTEDS
#jsonFileBUSTEDSMH  = snakemake.input.inputBUSTEDSMH
jsonFileBGM  = snakemake.input.inputBGM
jsonFileABSREL  = snakemake.input.inputABSREL
jsonFileSLAC  = snakemake.input.inputSLAC
jsonFileRELAX  = snakemake.input.inputRELAX
jsonFileCFEL  = snakemake.input.inputCFEL

jsons = {
         "FEL": jsonFileFEL,
         "FUBAR": jsonFileFUBAR,
         "MEME": jsonFileMEME,
         #"BUSTEDS": jsonFileBUSTEDS,
         #"BUSTEDSMH": jsonFileBUSTEDSMH,
         "BGM": jsonFileBGM,
         "ABSREL": jsonFileABSREL,
         "SLAC": jsonFileSLAC,
         "RELAX": jsonFileRELAX,
         "CFEL": jsonFileCFEL,
         }
         
#method = snakemake.params.method
output = snakemake.output.output
pvalueThreshold = 0.1
posteriorThreshold = 0.9

# Print to user

#print("# Summarizing", method, "JSON results")
#print("# Reading input from:", jsonFile)
print("# Saving output to:", output)

# Helper functions

def getJsonData(jsonFile: str):
    with open(jsonFile, "r") as fh:
        json_data = json.load(fh)
    #end with
    return json_data
#end method

def writeOutput(df, output):
    if not os.path.isfile(output):
        df.to_csv(output, index=False)
    else:
        df.to_csv(output, mode='a', index=False, header=False)
    #end if
#end method

# Main subroutine
# Main
def process(method, jsonFile, first=False):
    data = getJsonData(jsonFile)
    basename = os.path.basename(jsonFile)
    match method:
        case "FEL":
            FEL_data = data
            df = pd.DataFrame(FEL_data["MLE"]["content"]["0"],
                              columns=[x[0] for x in FEL_data["MLE"]["headers"]],
                              dtype = float)
            unadjusted_pvalues = df["p-value"].tolist()
            adjusted_pvalues =  statsmodels.stats.multitest.fdrcorrection(unadjusted_pvalues,
                                                                          alpha=0.10,
                                                                          method='indep',
                                                                          is_sorted=False)
            df["adjusted_p-value"] = adjusted_pvalues[1]
            df_results = df[df["adjusted_p-value"] <= pvalueThreshold]
            positive_sites = df_results[df_results["dN/dS MLE"] > 1.0]
            negative_sites = df_results[df_results["dN/dS MLE"] < 1.0]
            N = int(FEL_data["input"]["number of sites"])
            S = int(FEL_data["input"]["number of sequences"])
            results_dict  = {
                            "Filename": str(basename),
                            "Seqs": str(S),
                            "Sites": str(N),
                            "FEL[+]": str(positive_sites.shape[0]),
                            "FEL[-]": str(negative_sites.shape[0]),
                            }
            #df_output = pd.DataFrame.from_dict(results_dict, orient="index").T
            #writeOutput(df_output, output)
            
            return results_dict
        case "FUBAR":
            columns = data["MLE"]["headers"]
            headers = [x[0] for x in columns]
            headers.append("Z") # Placeholders
            headers.append("Y") # Placeholders
            df = pd.DataFrame(data["MLE"]["content"]["0"], columns=headers)
            positive_sites = df[df["Prob[alpha<beta]"] >= posteriorThreshold]
            negative_sites = df[df["Prob[alpha>beta]"] >= posteriorThreshold]
            results_dict  = {
                            "Filename": str(basename),
                            "FUBAR[+]": str(positive_sites.shape[0]),
                            "FUBAR[-]": str(negative_sites.shape[0]),
                            }
            return results_dict
            #df_output = pd.DataFrame.from_dict(results_dict, orient="index").T
            #writeOutput(df_output, output)
        case "MEME":
            MEME_data = data
            df = pd.DataFrame(MEME_data["MLE"]["content"]["0"], columns=[x[0] for x in MEME_data["MLE"]["headers"]], dtype = float)
            unadjusted_pvalues = df["p-value"].tolist()
            adjusted_pvalues =  statsmodels.stats.multitest.fdrcorrection(unadjusted_pvalues,
                                                                          alpha=0.10,
                                                                          method='indep',
                                                                          is_sorted=False)
            df["adjusted_p-value"] = adjusted_pvalues[1]
            df_results = df[df["adjusted_p-value"] <= pvalueThreshold]
            df_results = df[df["p-value"] <= pvalueThreshold]
            basename = os.path.basename(jsonFile)
            MEME_results = df_results.shape[0]
            results_dict  = {
                            "Filename": str(basename),
                            "MEME": str(df_results.shape[0])
                            }
            return results_dict
        case "BUSTEDS":
            # BUSTED[S]
            BUSTEDS_data = data
            BUSTEDS_pvalue = BUSTEDS_data["test results"]["p-value"]
            #BUSTED[S]+MH
            #print ("# Processing:", BUSTEDS_JSON)
            #BUSTEDSMH_data = get_JSON(BUSTEDSMH_JSON)
            #BUSTEDSMH_pvalue = BUSTEDSMH_data["test results"]["p-value"]
            # "BUSTED[S]": BUSTEDS_pvalue,
            results_dict  = {
                            "Filename": str(basename),
                            "BUSTED[S]": str(BUSTEDS_pvalue)
                            }
            return results_dict
        case "BUSTEDSMH":
            BUSTEDSMH_data = data
            BUSTEDSMH_pvalue = BUSTEDSMH_data["test results"]["p-value"]
            results_dict  = {
                            "Filename": str(basename),
                            "BUSTED[S]+MH": str(BUSTEDSMH_pvalue)
                            }
            return results_dict
            
        case "BGM":
            #posteriorThresholdBGM = 0.5
            #columns = data["MLE"]["headers"]
            #headers = [x[0] for x in columns]
            
            #print(headers)
            #headers2= []
            #for item in headers:
            #    item = item.replace('â€“', "-")
            #    headers2.append(item)
            #end for
            
            #df = pd.DataFrame(data["MLE"]["content"], columns=headers)
            
            #print(df)
            #df.to_csv("bgm_test.csv")
            #coevolving_sites_1 = df[df["P [Site 1 –> Site 2]"] >= posteriorThreshold]
            #coevolving_sites_2 = df[df["P [Site 2 -> Site 1]"] >= posteriorThreshold]
            
            #for col in df.columns:
            #    print([col], df[col])
            
            #coevolving_sites_3 = df[df["P [Site 1 <–> Site 2]"] >= posteriorThresholdBGM]
            #return coevolving_sites_3.shape[0]
            results_dict  = {
                            "Filename": str(basename),
                            "BGM": "0"
                            }
            return results_dict
        case "ABSREL":
            results_dict  = {
                            "Filename": str(basename),
                            "ABSREL": "0"
                            }
            return results_dict
        case "SLAC":
            results_dict  = {
                            "Filename": str(basename),
                            "SLAC": "0"
                            }
            return results_dict
        case "RELAX":
            results_dict  = {
                            "Filename": str(basename),
                            "RELAX": "0"
                            }
            return results_dict
        case "CFEL":
            results_dict  = {
                            "Filename": str(basename),
                            "CFEL": "0"
                            }
            return results_dict
        case _:
            pass
        #end case
    #end match
#end method

results_dict = {}

for _ in jsons.keys():
    print("# Processing...", _, jsons[_])
    if results_dict == {}:
        results_dict = process(_, jsons[_], True)
    else:
        # d1.update(d2)
        results_dict.update(process(_, jsons[_]))
    #end if
#end for

df_output = pd.DataFrame.from_dict(results_dict, orient="index").T
writeOutput(df_output, output)


# End of file
