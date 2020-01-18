import pandas as pd
import numpy as np

def format_genes(gene_coord_df, chr_sizes_df, out_dir, windowsize_kb, print_log_files=True):
	""" 
	Adds a fixed-size window to each protein-coding gene in the human genome.

	Input
		gene_coord_df: Locations of all human genes in hg19
		chr_sizes_df: Sizes of hg19 human autosomes and sex chromosomes
		out_dir: Where the output BED file will be saved
	Output

	"""
	# A dictionary where keys are chromosomes, values are the lengths in BP
	chr_sizes_dict = chr_sizes_df.to_dict()[1]

	# Gets the window size in bases
	windowsize = windowsize_kb * 1000

	# Keep only protein-coding genes, and genes on autosomes or sex chromosomes 
	# (Remove sex chromosomes later in pipeline as LDSC only uses autosomes for now)
	gene_coord_df = gene_coord_df[gene_coord_df["gene_biotype"]=='protein_coding']
	gene_coord_df = gene_coord_df[gene_coord_df['CHR'].isin(chr_sizes_dict.keys())]
	gene_coord_df = gene_coord_df.dropna().sort_values(by=['CHR','START'])


	gene_coord_df = gene_coord_df.drop('gene_biotype', 1)
	gene_coord_df['START'] = np.maximum(0, gene_coord_df['START'] - windowsize)
	gene_coord_df['END'] = np.minimum(gene_coord_df['END'] + windowsize, gene_coord_df['CHR'].map(chr_sizes_dict))

	# Save
	for chrom_num in chr_sizes_dict.keys():
		chrom_df = gene_coord_df[gene_coord_df['CHR'] == chrom_num]
		chrom_df = chrom_df[['CHR','START','END','GENE']]
		chrom_df.to_csv('{out_dir}/genes_plus_{window}kb.{chr}.bed'.format(out_dir = out_dir,
																			window = windowsize_kb,
																			chr = chrom_num), sep='\t', index=False, header = False)


###################################### MAIN ######################################

snake_log_obj = snakemake.log # class(snakemake.log) = 'snakemake.io.Log
sys.stdout = open(str(snake_log_obj), "w") # could not find better ways than calling str(snake_log_obj) to get the log filepath


windowsize_kb = snakemake.params['windowsize_kb']
bed_out_dir = snakemake.params['bed_out_dir']


gene_coords = snakemake.input['gene_coords']
chr_sizes = snakemake.input['chr_sizes']

gene_coords_df = pd.read_csv(gene_coords, delim_whitespace = True, index_col=None)
chr_sizes_df = pd.read_csv(chr_sizes, delim_whitespace = True, index_col=0, header=None)
format_genes(gene_coords_df, chr_sizes_df, bed_out_dir, windowsize_kb)
# multi_gene_sets_to_dict_of_beds(df_multi_gene_set_human, df_gene_coords, windowsize, bed_out_dir + '/tmp', bed_out_dir, out_prefix)