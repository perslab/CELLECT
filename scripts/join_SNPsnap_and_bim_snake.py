import pandas as pd
import numpy as np
import argparse
import os
from time import gmtime, strftime
import multiprocessing
import csv


def join_SNPsnap_bim_SNPs_to_input_bim(input_bim_path, output_dir, SNPsnap_file, chromosomes):
	'''
	Left joins the input_bim file with the SNPsnap file then splits the resulting 
	joined dataframe by chromosome and saves them each as csv files.

	Args:
		input_bim_path (str): A path to the input bim file containing the SNPs 
		we wish to analyse.
		SNPsnap_path (str): A path to the SNPsnap file containing the SNPs and 
		other information.

	Returns:
		None
	'''
	# Reading in both bim files as panda dataframes

	SNPsnap_SNPs = pd.read_csv(SNPsnap_file, sep='\t')
	print('There are {} SNPs in the SNPsnap file'.format(len(SNPsnap_SNPs)))

	full_bim_path = input_bim_path + '.CHR_1_22.bim'
	input_SNPs = pd.read_csv(full_bim_path, sep='\t', header=None,
	 names=['Chromosome',
	  'rsID', 
	  'Position in morgans', 
	  'BP coordinate', 
	  'Allele 1', 
	  'Allele 2'])
	print('There are {} SNPs in the input bim file'.format(len(input_SNPs)))
	input_bim_filename = os.path.basename(input_bim_path)

	# Making a new column to join SNPs on as rsIDs don't match
	input_SNPs['snpID'] = input_SNPs.loc[:,'Chromosome'].astype(str) + ':' + input_SNPs.loc[:,'BP coordinate'].astype(str)


	# Joining the SNPsnap snps onto the input SNPs
	print('Left joining input bim file with the SNPsnap file.')
	left_joined_SNPs = input_SNPs.join(SNPsnap_SNPs.set_index('snpID'),
		on='snpID', lsuffix='_in', rsuffix='_snap', how='left')
	print('{} SNPs were lost in the joining process'.format(np.sum(left_joined_SNPs.loc[:,'rsID_snap'].isna())))

	# Dropping columns preparing for saving the final output
	left_joined_SNPs = left_joined_SNPs.loc[:,['snpID', 'ID_genes_in_matched_locus', 'Chromosome', 'rsID_in', 'rsID_snap', 'BP coordinate']]
	left_joined_SNPs.sort_values(by='snpID', inplace=True, ascending=True)

	for chr_num in chromosomes: 
		cur_chr = left_joined_SNPs[left_joined_SNPs.loc[:,'Chromosome']==chr_num]
		cur_chr_out_path = '{output_dir}/SNPs_with_genes.{bim_name}.{chromosome}.txt'.format(output_dir = output_dir,
																						bim_name = input_bim_filename,
																						chromosome = chr_num)
		cur_chr.to_csv(cur_chr_out_path, sep='\t', index=False)



	print('SNPs with gene files successfully produced for chromosomes {chromosome_numbers} in {bim_name}'.format(bim_name = input_bim_filename,
																												chromosome_numbers = ', '.join(map(str, chromosomes))))
	return(None)


###################################### MAIN ######################################

bimfile_path = snakemake.config['LDSC_BFILE_PATH']
SNPsnap_file = snakemake.config['LDSC_SNPSNAP_FILE']

out_dir = snakemake.params['out_dir']
list_chromosomes = snakemake.params['chromosomes']

join_SNPsnap_bim_SNPs_to_input_bim(bimfile_path, out_dir, SNPsnap_file, list_chromosomes)