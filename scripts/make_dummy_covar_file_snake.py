import pandas as pd
import numpy as np



### Enable logging
snake_log_obj = snakemake.log # class(snakemake.log) = 'snakemake.io.Log
sys.stdout = open(str(snake_log_obj), "w") # could not find better ways than calling str(snake_log_obj) to get the log filepath

### Get parameters from snakemake
GENELOC_FILE = snakemake.params['geneloc_file']
DUMMY_COVAR_FILE = snakemake.params['dummy_covar_file']

print("Reading gene locations from the file " + str(GENELOC_FILE))
all_magma_genes = pd.read_csv(GENELOC_FILE, sep='\t', header = None)

print("Composing a dummy covariance file...")
all_magma_genes = all_magma_genes.rename(columns={ all_magma_genes.columns[0]: "entrez_id" })
all_magma_genes = pd.DataFrame(all_magma_genes.iloc[:, 0])
all_magma_genes['dummy_gene_covar'] = np.random.normal(size = len(all_magma_genes))

all_magma_genes.to_csv(DUMMY_COVAR_FILE, sep='\t', index = None)

print("The dummy covar file " + str(DUMMY_COVAR_FILE) + " is created. The top 10 positions are: ")
print(all_magma_genes.head(10))

