import pandas as pd
import numpy as np

### File with Entrez - Ensembl id pairs (a global variable)
mapping_file = snakemake.params['mapping_file']                      


### The mapping Function:  ENTREZ --> ENSEMBL
def add_ensembl_ids_from_entrez(df, colname_geneids_from = 'entrez', colname_geneids_to = 'ensembl_gene_id'):
	### INPUT df: a data frame with the column 'colname_geneids_from' with entrez gene ids.
	### OUTPUT df:
	# returns a data frame with ensembl gene ids added to the column 'colname_geneids_to'.
	# Genes that did not map have NA values in the 'colname_geneids_to' column.
	# If there are duplicated gene IDs in 'colname_geneids_to', then all but the first of the duplicated elements will be marked as NA.
	# 'mapping_file' containing the corresponding  pais of Ensembl - Entrez gene ids is used for mapping.        

	df_mapping = pd.read_csv(mapping_file, sep = '\t')
	merged_df = pd.merge(df, df_mapping, left_on = 'GENE', right_on = 'entrezgene', how = 'inner')
	merged_df = merged_df.drop_duplicates('entrezgene')
	merged_df = merged_df.drop(columns = ['entrezgene'])
	print("Number of genes mapped: " + str(len(merged_df)))
	unmapped_genes = df[df[colname_geneids_from].isin(df_mapping['entrezgene']) == False]
	print("Number of genes not mapped: " + str(len(unmapped_genes)))
	merged_df.loc[merged_df.duplicated('ensembl_gene_id'), 'ensembl_gene_id'] = np.nan # find elements with smaller subscripts, mark them as NA (duplicates!)
	dups = merged_df[merged_df['ensembl_gene_id'].isnull()]
	print("Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): " + str(len(dups)))
	# concat mapped and unmapped genes into one dataframe
	res = pd.concat([merged_df, unmapped_genes], sort = True)
	res = res.rename(columns = {"ensembl_gene_id": "gene"})
	print("Total genes mapped (non NA genes): " + str(len(merged_df) - len(dups)))
	print("Returning dataframe with the column '" + colname_geneids_to +  "' added where all gene identifiers are unique. Unmapped genes have NA values.")
	return res





########################################## MAIN ##########################################

### Enable logging
snake_log_obj = snakemake.log                  # class(snakemake.log) = 'snakemake.io.Log
sys.stdout = open(str(snake_log_obj), "w")     # could not find better ways than calling str(snake_log_obj) to get the log filepath

### Get the other parameters from snakemake
gwas_name = snakemake.params['gwas_name']                           # GWAS name
base_output_dir = snakemake.params['base_output_dir']               # base directory for all the outputs

### Open MAGMA ZSTATs file
df_magma = pd.read_csv(base_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".resid_correct_all.gsa.genes.out", comment = '#', sep = '\s+', header = 0)

print("Mapping human ENTREZ gene IDs to human ENSEMBL gene IDs for " + gwas_name + ": \n")

### Mapping: Entrez gene IDs --> Ensembl gene IDs
df_magma = add_ensembl_ids_from_entrez(df_magma, colname_geneids_from = 'GENE', colname_geneids_to = 'gene')

### Exclude unmapped genes
print("Excluding unmapped genes...")
df_magma = df_magma.loc[:, df_magma.columns != 'GENE']
df_magma = df_magma.dropna(subset = ['gene'])

### Save the result for the given GWAS
print("Saving the results...")
res_fullname = base_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".resid_correct_all_ens.gsa.genes.out"
df_magma.to_csv(res_fullname, sep = '\t', index = False)
print("The results are saved as " + res_fullname)
