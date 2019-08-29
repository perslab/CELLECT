import os
import sys
import re
import pandas as pd
import pdb

# import glob
# import argparse
# import subprocess
# import multiprocessing

def split_M_files(file_ldscore, chromosome, prefix_genomic_annot, list_annotations):
	### Both file_M and file_M_5_50 contains a single line with the same number of fields as the number of annotations in the ldscore file
	# nn_lira_sema.COMBINED_ANNOT.6.l2.ldscore.gz
	# nn_lira_sema.COMBINED_ANNOT.6.l2.M
	# nn_lira_sema.COMBINED_ANNOT.6.l2.M_5_50
	
	print("CHR={} | Writing .M and .M_5_50 files".format(chromosome))
	file_ldscore_base = re.sub(r"\.l2\.ldscore\.gz$", "", file_ldscore) # .e.g /scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.COMBINED_ANNOT.9.
	file_M = "{}.l2.M".format(file_ldscore_base)
	file_M_5_50 = "{}.l2.M_5_50".format(file_ldscore_base)
	with open(file_M, "r") as fh_M, open(file_M_5_50, "r") as fh_M_5_50:
		list_M =  fh_M.readline().rstrip().split()
		list_M_5_50 = fh_M_5_50.readline().rstrip().split()
	assert(len(list_M) == len(list_M_5_50) == len(list_annotations))
	for i in range(len(list_annotations)):
		annotation = list_annotations[i]
		M = list_M[i]
		M_5_50 = list_M_5_50[i]
		file_out_M = "{}/{}__{}.{}.l2.M".format(out_dir, prefix_genomic_annot, annotation, chromosome)
		file_out_M_5_50 = "{}/{}__{}.{}.l2.M_5_50".format(out_dir, prefix_genomic_annot, annotation, chromosome)
		with open(file_out_M, "w") as fh_out_M, open(file_out_M_5_50, "w") as fh_out_M_5_50:
			fh_out_M.write(M)
			fh_out_M_5_50.write(M_5_50)
	print("CHR={} | DONE writing .M and .M_5_50 files".format(chromosome))


def split_annot_file(file_ldscore, chromosome, prefix_genomic_annot, list_annotations):
	""" Split annot file 
	Because list_annotations is looped over, this function should work for both thin-annot and standard annot.
	"""
	print("CHR={} | START splitting annot file".format(chromosome))
	file_ldscore_base = re.sub(r"\.l2\.ldscore\.gz$", "", file_ldscore) # .e.g /scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.COMBINED_ANNOT.9.
	file_annot = "{}.annot.gz".format(file_ldscore_base)
	df_annot = pd.read_csv(file_annot, sep="\t") # no index, compression detected automatically
	print("CHR={} | Read  annot file: {}".format(chromosome, file_annot))
	for i in range(len(list_annotations)):
		if i % 25 == 0:
			print("CHR={} #{}/#{}| Writing annot files".format(chromosome, i, len(list_annotations)))
		annotation = list_annotations[i]
		file_out_annot = "{}/{}__{}.{}.annot.gz".format(out_dir, prefix_genomic_annot, annotation, chromosome)
		df_annot[[annotation]].to_csv(file_out_annot, sep="\t", index=False, compression="gzip") # df_annot[[annotation]] to return a DataFrame with a header instead of Series with no header (df_annot[annotation] returns this).
	print("CHR={} | DONE splitting annot file".format(chromosome))

def split_ldscore_file_per_annotation(file_ldscore):
	print("Processing file_ldscore {}".format(file_ldscore))
	m = re.search(r"(.*)\.COMBINED_ANNOT\.(\d{1,2})\.l2.ldscore.gz$", os.path.basename(file_ldscore)) # file_ldscore e.g. "/scratch/sc-ldsc/nn_lira_sema/nn_lira_sema.COMBINED_ANNOT.9.l2.ldscore.gz"
	prefix_genomic_annot = m.groups()[0] # e.g nn_lira_sema
	chromosome = m.groups()[1]
	df = pd.read_csv(file_ldscore, sep="\t") # no index
	# CHR     SNP     BP      blueL2  cornsilkL2 ...
	annotations_header = df.columns[3:].tolist() # skip the first three columns (CHR, SNP, BP)
	annotations_clean = [re.sub(r"L2$", "", x) for x in annotations_header] # remove trailing L2 in name. Do not use x.rstrip("L2")
	split_M_files(file_ldscore, chromosome, prefix_genomic_annot, annotations_clean)
	split_annot_file(file_ldscore, chromosome, prefix_genomic_annot, annotations_clean)
	for counter, annotation in enumerate(annotations_header):
		if counter % 25 == 0:
			print("CHR={} #{}/#{}| Writing ldscore files".format(chromosome, counter, len(annotations_header)))
		annotation_clean = re.sub(r"L2$", "", annotation)
		file_out_ldscore = "{}/{}__{}.{}.l2.ldscore.gz".format(out_dir, prefix_genomic_annot, annotation_clean, chromosome)
		df[["CHR", "SNP", "BP", annotation]].to_csv(file_out_ldscore, sep="\t", index=False, compression="gzip") # no index



###################################### MAIN ######################################

# snake_log_obj = snakemake.log # class(snakemake.log) = 'snakemake.io.Log
# sys.stdout = open(str(snake_log_obj), "w") # could not find better ways than calling str(snake_log_obj) to get the log filepath


file_ldscore = snakemake.input["ldscore"] # e.g. <SOMEPATH>/nn_lira_sema.COMBINED_ANNOT.22.l2.ldscore.gz
chromosome = snakemake.params['chromosome']
dataset_prefix = snakemake.params['run_prefix']

out_dir = os.path.join(snakemake.params['out_dir'], "per_annotation")
split_ldscore_file_per_annotation(file_ldscore)









