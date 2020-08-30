import pandas as pd
import numpy as np
import os
import scipy.stats as st

def get_pvals(df_magma):
    '''
    usage: convert Z stats to p-values
    '''
    df_magma["P"] = 1-st.norm.cdf(df_magma["ZSTAT"]) # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html
    # subtracting an array from 1 gives an array
    return df_magma


########################################################## MAIN ###############################################################

### Get the parameters from snakemake
gwas_name = snakemake.params['gwas_name']				   # GWAS name
#base_output_dir = snakemake.params['base_output_dir']			   # base directory for all the outputs

base_output_dir = snakemake.params['base_output_dir']

print("Computing corrected p-values for GWAS '" + gwas_name + "': ")

### Load df_magma
df_magma = pd.read_csv(base_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".resid_correct_all.gsa.genes.out", sep= '\s+', header = 1)

df_out = get_pvals(df_magma)

### Save the results for the given specificity_id and GWAS
outname = gwas_name + "/" + gwas_name + ".resid_correct_all.gsa.genes.pvals.out"
outdir = base_output_dir + "/precomputation/"

fullname = os.path.join(outdir, outname)

print("Saving the results...")
df_out.to_csv(fullname, sep= '\t', index = False)
print("The results are saved as " + fullname)
