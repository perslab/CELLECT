import pandas as pd
import numpy as np


###################################### FUNCTIONS ######################################

def make_all_genes(specificity_matrix_file, specificity_matrix_name, out_dir):
	specificity_df = pd.read_csv(specificity_matrix_file, index_col='gene')

	all_genes = pd.DataFrame(index = specificity_df.index)
	all_genes['all_genes_in_dataset'] = 1

	all_genes_path = '{out_dir}/all_genes.{sm_name}.csv'.format(
						out_dir = out_dir,
						sm_name = specificity_matrix_name)
	all_genes.to_csv(all_genes_path)


###################################### MAIN ######################################

output_dir = snakemake.params['out_dir']
specificity_matrix_name = snakemake.params['specificity_matrix_name']
specificity_matrix_file = snakemake.params['specificity_matrix_file']

make_all_genes(specificity_matrix_file, specificity_matrix_name, output_dir)