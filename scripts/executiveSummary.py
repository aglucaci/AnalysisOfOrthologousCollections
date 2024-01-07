"""

"""

import pandas as pd
import numpy as np
import os
import json
import glob
from os.path import exists
import sys

# Declares


jsonFile = snakemake.input.FEL_JSON
#outputPNG = snakemake.output.output_png
#outputCSV = snakemake.output.output_csv
#outputFigureLegend = snakemake.output.figure_legend
outputCumulativeResults = snakemake.output.cumulative_results

pvalueThreshold = 0.1

method = snakemake.params.METHOD

print("# Summarizing", method, "JSON results")

data = {}

# Helper functions

def helper():
    pass
#end method

# Main 

data[0] = {"TestCol": "TestValue", 
           "Filename": jsonFile,
           "Method": method,
           "CodonSite": "1",
           "p-value": 0,
           "dNdS":0}

df = pd.DataFrame(data)

# Save the DataFrame to a CSV file
df.to_csv(outputCumulativeResults, index=True)

print(f"DataFrame saved as '{outputCumulativeResults}'")

# End of file
