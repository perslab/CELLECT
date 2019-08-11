import pandas as pd
import numpy as np

out_dir = snakemake.params['out_dir']
specificity_matrix_name = snakemake.params['specificity_matrix_name']
specificity_matrix_file = snakemake.params['specificity_matrix_file']


specificity_df = pd.read_csv(specificity_matrix_file)

# Convert the specificty dataframe into the long multigeneset format
# required as input for LD score regression
multi_geneset = pd.melt(specificity_df,
							 id_vars=['gene'],
							 var_name='annotation',
							 value_name='specificity')
multi_geneset = multi_geneset[['annotation','gene','specificity']]
multi_geneset = multi_geneset.loc[multi_geneset.specificity>0]


# Make a multigeneset from all genes in the dataset to serve as a background
multi_geneset_all_genes = pd.DataFrame(data={'annotation':'all_genes_in_dataset',
											 "gene":np.unique(specificity_df.gene),
											 'specificity':1})

# Save both multigeneset dataframes as tab separated files
mgs_path = '{out_dir}/multi_geneset.{sm_name}.txt'.format(
			out_dir = output_dir,
			sm_name = specificity_matrix_name)
multi_geneset.to_csv(mgs_path,
					 header=None, index=False,sep='\t')

mgs_all_genes_path = '{out_dir}/multi_geneset.{sm_name}.all_genes.txt'.format(
					  out_dir = output_dir,
					  sm_name = specificity_matrix_name)
multi_geneset_all_genes.to_csv(mgs_all_genes_path,
							   header=None, index=False,sep='\t')