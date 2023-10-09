#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import os
import json
import altair as alt
import altair_saver
from altair_saver import save
import argparse

# In[3]:


#JSON_FILE = snakemake.input.input
#OUTPUT_PNG = snakemake.output.output_png
#OUTPUT_CSV = snakemake.output.output_csv
#output_figure_legend = snakemake.output.figure_legend

parser = argparse.ArgumentParser(
                    prog='FEL Summary',
                    description='Summarizes FEL results',
                    epilog='This script processes FEL results, no statistical adjustment')

parser.add_argument('--json-file', required=True)           
parser.add_argument('--output-png', required=True)      
parser.add_argument('--output-csv', required=True) 
parser.add_argument('--output-figurelegend', required=True) 

args = parser.parse_args()
#print(args.filename, args.count, args.verbose)

JSON_FILE = args.json_file
OUTPUT_PNG = args.output_png
OUTPUT_CSV = args.output_csv
output_figure_legend = args.output_figurelegend

pvalueThreshold = 0.1


# In[4]:


def getFELData(json_file):
    with open(json_file, "r") as in_d:
        json_data = json.load(in_d)
    return json_data["MLE"]["content"]["0"]
#end method

def getFELHeaders(json_file):
    with open(json_file, "r") as in_d:
        json_data = json.load(in_d)
    return json_data["MLE"]["headers"]
#end method


# In[5]:


columns = getFELHeaders(JSON_FILE)
headers = [x[0] for x in columns]


# In[21]:


data = getFELData(JSON_FILE)


# ### Selected Sites -- Tables

# In[22]:


df = pd.DataFrame(getFELData(JSON_FILE), columns=headers, dtype = float)
df["omega"] = df["beta"] / df["alpha"]
df.index += 1
df["Site"] = df.index

# Saving CSV
df.to_csv(OUTPUT_CSV)


# In[ ]:


df_results = df[df["p-value"] <= pvalueThreshold]


# In[24]:


positive_sites = df_results[df_results["omega"] > 1.0]
positive_sites = positive_sites.reset_index()
positive_sites.index += 1
positive_sites.drop('index', axis=1, inplace=True)
#positive_sites


# In[25]:


negative_sites = df_results[df_results["omega"] < 1.0]
negative_sites = negative_sites.reset_index()
negative_sites.index += 1
negative_sites.drop('index', axis=1, inplace=True)
#negative_sites


# ## Visualizations

# In[ ]:


source = df

line = alt.Chart(source).mark_line().encode(
    x='Site',
    y = alt.Y('dN/dS MLE')
).properties(
    width=800,
    height=600)


band = alt.Chart(source).mark_area(opacity=0.5).encode(x='Site',
                                                       y='dN/dS LB', 
                                                       y2='dN/dS UB')

chart = line+band

#save(chart, OUTPUT_PNG)

print("# Saving chart:", [OUTPUT_PNG])

save(chart, OUTPUT_PNG)

#line + band
# altair saver


# ## Figure legend.

# In[2]:


## Summary

a = len(df["omega"])
b = len(negative_sites["omega"])
d = len(positive_sites["omega"])
c = round((b/a) * 100, 2)

pct_neg = c
pct_pos =  round((d/a) * 100, 2)


with open(output_figure_legend, "w") as fh:
    print("FEL analysis of your gene of interest found " + str(b) + " of " + str(a) + 
          " (" + str(pct_neg) + "%) sites to be statistically significant (p-value ≤ " + str(pvalueThreshold) + 
          ") for pervasive negative/purifying selection. In addition, we observe evidence that " + str(d) + " (" + str(pct_pos) +
          "%) sites are operating under a positive/adaptive selection regime.", file=fh)

