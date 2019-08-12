import multiprocessing
import os
import sys
import subprocess

import itertools
import functools

import pandas as pd
import numpy as np
import argparse

import pybedtools # this script was developed using pybedtools v. 0.8.0
# ^ some/all of pybedtools requires 'bedtools' to be available on your PATH /tools/bedtools/2.27.1/bin/bedtools
# BedTool(..).sort(): "sortBed" does not appear to be installed or on the path, so this method is disabled. Please install a more recent version of BEDTools and re-import to use this method.
# Install via conda install --channel conda-forge --channel bioconda pybedtools bedtools htslib
# [did not work for me, so I used conda install ... instead] Use the below lines to add a system installation of bedtools and tabix when running script within anaconda:
# os.environ["PATH"] += os.pathsep + "/tools/bedtools/2.27.1/bin/"
# os.environ["PATH"] += os.pathsep + "/tools/htslib/1.6/bin/"

import psutil





###################################### DOCUMENTATION ######################################

### General remarks:
# - whitespace or forward slash are not allowed in annotation_name column in file_multi_gene_set.
#   This is because LDSC .annot files are read as *whitespace delimted* by the ldsc.py program, so annotation_name with whitespace in the name will make the .l2.ldscore.gz header wrong.
# - If a log.<out_prefix>.multi_geneset.txt file exists, the file_multi_gene_set/df_multi_gene_set must contain the same annotations (or a subset of them) with the same annotation values.

### For continuous annotations:
# - When a variant is spanned by multiple genes with the XXX kb window, we assign the maximum annotation_value.

### For binary annotations:
# - The SNPs within the genomic regions spanned by the genes within a given annotation gets the annotation value 1. All other SNPs get the annotation value 0

### Resource usage:
# - This script is memory intensive for inputs with many annotations. See notes on "Determine number of parallel processes to run"
# - This script is I/O intensive because pybedtools with write many temporary files during Bedtools operations.


###################################### Functions ######################################

def make_annot_file_per_chromosome(chromosome, dict_of_beds, out_dir, out_prefix, annot_per_geneset, bimfile, all_genes):
	""" 
	Input
		chromosome: integer (1..22)
	
	*OBS* this function RELIES on MANY GLOBAL scope VARIABLES
	"""
	# TODO: parse variables to function
	
	### make annot file
	print('making annot files for chromosome {}'.format(chromosome))
	df_bim = pd.read_csv(bimfile, delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
	# (Pdb) df_bim.head()
	#    CHR          SNP        CM       BP
	# 0   21  rs146134162 -0.908263  9412099
	# 1   21  rs578050168 -0.908090  9412377
	# 2   21  rs527616997 -0.907297  9413645
	# 3   21  rs544748596 -0.906578  9414796
	# 4   21  rs528236937 -0.906500  9414921
	# iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
	# ^ Python3 (but not Python2.7) gives the following error when calling "bimbed = BedTool(iter_bim)" in make_annot_file_per_chromosome()
	# /tools/anaconda/3-4.4.0/lib/python3.6/site-packages/pybedtools/cbedtools.pyx in pybedtools.cbedtools.IntervalIterator.__next__()
	# /tools/anaconda/3-4.4.0/lib/python3.6/site-packages/pybedtools/cbedtools.pyx in pybedtools.cbedtools.create_interval_from_list()
	# /tools/anaconda/3-4.4.0/lib/python3.6/site-packages/pybedtools/cbedtools.pyx in pybedtools.cbedtools.isdigit()
	# AttributeError: 'numpy.int64' object has no attribute 'isdigit'
	# SOLUTION: convert everything to strings --> ['chr'+str(x1), str(x2), str(x2)]
	print(df_bim.head)
	print(bimfile)
	iter_bim = [['chr'+str(x1), str(x2), str(x2)] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
	bimbed = pybedtools.BedTool(iter_bim)
	counter = 1 # just to print status message
	list_df_annot = []
	for name_annotation in sorted(dict_of_beds): # we sort to make output more consistent.
		print("CHR={} | annotation={}, #{}/#{}".format(chromosome, name_annotation, counter, len(dict_of_beds)))
		bed_for_annot = dict_of_beds[name_annotation] # get bed
		# (Pdb)len(bed_for_annot)
		# 66
		#  (Pdb) bed_for_annot.head()
		# chr1    51619935        52185000
		#  chr1   70410488        70871303
		#  chr1   85584164        86243933
		#  chr1   202948059       203355877
		#  chr10  43851792        44270066
		#  chr10  75681524        76110821
		#  chr10  76769912        77191206
		#  chr10  120663598       121138345
		#  chr11  118030300       118469926
		#  chr12  21454715        21871342
		
		annotbed = bimbed.intersect(bed_for_annot, wb=True) # PT NOTE: this finds SNPs in bim file that OVERLAP with the annotation bed (gene)
		# chr22  24008141    24008141    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008403    24008403    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008409    24008409    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008465    24008465    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008495    24008495    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008497    24008497    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008503    24008503    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008699    24008699    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		# chr22  24008773    24008773    chr22   24008021    24210630    ENSG00000250479 0.03038823367
		
		# annotbed = bed_for_annot.intersect(bimbed, wb=True) # PT NOTE: this finds the positions/intervals in annotation bed (gene) that OVERLAP with the bim file. Only the part of the record intersections occurred
			# *IMPORTANT*: bimbed.intersect(bed_for_annot) and bed_for_annot.intersect(bimbed) DOES NOT return the same positions. However, they do return the same number of 'intersected features'. That is, the returned BedTool object as the same length.
			# bed_for_annot.intersect(bimbed) returns features that span two bp (e.g. start=24008140, end=24008142), whereas bimbed.intersect(bed_for_annot) returns features that span a single bp (start=24008141, end=24008141)
			# use bed_for_annot.intersect(bimbed, wb=True) to understand this behavior better.
		# chr22  24008140    24008142    ENSG00000250479 0.03038823367   chr22   24008141    24008141
		# chr22  24008402    24008404    ENSG00000250479 0.03038823367   chr22   24008403    24008403
		# chr22  24008408    24008410    ENSG00000250479 0.03038823367   chr22   24008409    24008409
		# chr22  24008464    24008466    ENSG00000250479 0.03038823367   chr22   24008465    24008465
		# chr22  24008494    24008496    ENSG00000250479 0.03038823367   chr22   24008495    24008495
		# chr22  24008496    24008498    ENSG00000250479 0.03038823367   chr22   24008497    24008497
		# chr22  24008502    24008504    ENSG00000250479 0.03038823367   chr22   24008503    24008503
		# chr22  24008698    24008700    ENSG00000250479 0.03038823367   chr22   24008699    24008699
		# chr22  24008772    24008774    ENSG00000250479 0.03038823367   chr22   24008773    24008773
		
		### DOCS .intersect()
		# the intervals reported are NOT the original gene intervals, but rather a refined interval reflecting solely the portion of each original gene interval that overlapped with the SNPs
		# The -wa (write A) and -wb (write B) options allow one to see the original records from the A and B files that overlapped. As such, instead of not only showing you where the intersections occurred, it shows you what intersected.
		# SEE MORE HERE: http://quinlanlab.org/tutorials/bedtools/bedtools.html
		# SEE https://daler.github.io/pybedtools/intersections.html
		# SEE https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.intersect.html#pybedtools.bedtool.BedTool.intersect        
		
		### Extract data from annotbed before deleting annotbed.fn
		### These iterations REQUIRE the tmp bed file (annotbed.fn) to exist. SEE cbedtools.IntervalFile or .cbedtools.IntervalIterator,
		### Since the iterables in annotbed (x.start or .fields[7]) are all strings (or integers?), they are immutable objects and python will insert the value of the object/string, not a reference, in the list comprehension
		### CONCLUSION: pass-by-value of immutable objects means that we can SAFELY DELETE annotbed.fn after this data.
		### REF Parameter Passing for Mutable & Immutable Objects: https://medium.com/@tyastropheus/tricky-python-ii-parameter-passing-for-mutable-immutable-objects-10e968cbda35
		bp = [x.start for x in annotbed] # PT NOTE: make list of all bp positions for the overlapping SNPs | All features, no matter what the file type, have chrom, start, stop, name, score, and strand attributes.
		annotation_value = [x.fields[7] for x in annotbed] # returns list of strings. Extract the 'score' column. This is column 7 in the 0-based column indexing. *OBS*: x.fields[7] is a string.
		
		### pybedtools cleanup V1: deletes all pybedtools session files [does not work - see below]
		### KEEP THIS AS A WIKI/EXPLANATION
		### REF 1 Pybedtools Design principles: https://daler.github.io/pybedtools/topical-design-principles.html
		### REF 2 https://daler.github.io/pybedtools/autodocs/pybedtools.helpers.cleanup.html#pybedtools.helpers.cleanup
		# Using BedTool instances typically has the side effect of creating temporary files on disk: every BedTools operation results in a new temporary file.
		# Temporary files may be created in order to run BEDTools programs, so BedTool objects must always point to a file on disk. 
		# Temporary files are stored in /tmp by default, and have the form /tmp/pybedtools.*.tmp.
		# By default, at exit all temp files created during the session will be deleted. 
		# However, if Python does not exit cleanly (e.g., from a bug in client code), then the temp files will not be deleted.
		# print("CHR={} | annotation={}, #{}/#{}. Doing pybedtools cleanup...".format(chromosome, name_annotation, counter, len(dict_of_beds)))
		# pybedtools.cleanup(verbose=True) # force deletion of all temp files from the current session.
		#  ---> YOU CANNOT CLEAN UP FILES at this point because it REMOVES ALL tmp files in dict_of_beds.
		#  ---> e.g. you get the execption: pybedtools.cbedtools.BedToolsFileError: /tmp/pybedtools.izu3ifzg.tmp does not exist

		### pybedtools cleanup V2: deletes current annotbed (specific to a chromosome and annotation)
		# we need to cleanup files because a lot of tmp bed files are written to pybedtools.get_tempdir().
		# tmp bed files can take up to 200 MB per file. The file size is dependent on the number of genes in the annotation. 
		# So inputs with "raw SEMs" annotations where each annotation contains all genes in the dataset (all genes have a non-zero SEM value) will generate large tmp bed files.
		# >>900 GB storage is used if running 2-4 parallel processes of make_annot_file_per_chromosome() and ~1500 annotations
		# Summary of storage use for this function if not doing forced clean-up : N_files = N_annotations * N_parallel_procs. e.g. 1500 annotations * 4 proc * 0.2 GB per file = 1200 GB
		# By default, tmp bed files would only cleaned up after completing this function. 
		# OUR SOLUTION: after doing the Bedtools intersect opertation, we no longer need the tmp file (the annotbed object lives in python memory). Force removal of the tmp bed file specific to a chromosome and annotation
		# NOTE: deleting a tmp file will not cause any problems later on for pybedtools automatic cleanup. I check the source code.
		os.remove(annotbed.fn)

		### Create data frame
		df_annot_overlap_bp = pd.DataFrame({'BP': bp, name_annotation:annotation_value}) # FINUCANE ORIG: df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
		#             BP  blue
		# 0     34605531     1
		# 1     34605604     1
		# 2     34605644     1
		# 3     34605778     1
		# 4     34606634     1
		# 5     34606840     1
		# 6     34607223     1
		df_annot = pd.merge(df_bim, df_annot_overlap_bp, how='left', on='BP') # *IMPORTANT*: how='left' --> resulting data frame will include ALL snps from the bim file.
		# ^ how="left": use only keys from left frame, PRESERVE KEY ORDER
		# (Pdb) df_annot.head()
		#    CHR          SNP        CM       BP  blue
		# 0   21  rs146134162 -0.908263  9412099   NaN
		# 1   21  rs578050168 -0.908090  9412377   NaN
		# 2   21  rs527616997 -0.907297  9413645   NaN
		# 3   21  rs544748596 -0.906578  9414796   NaN
		# 4   21  rs528236937 -0.906500  9414921   NaN
		df_annot = df_annot[[name_annotation]] # get rid of all columns but the name_annotation. Important: return 1 column data frame (and not series, which would loose the column name)
			# df[[name_annotation]] or df.loc[:, [name_annotation]] --> returns dataframe
			# df[name_annotation] or df.loc[:, name_annotation] --> returns series
		df_annot.fillna(0, inplace=True) # SNPs not in df_annot_overlap_bp will have NA values in name_annotation
		# Do data type conversion AFTER .fillna() to avoid problems with NA (float) that cannot be converted to int.
		df_annot[name_annotation] = df_annot[name_annotation].astype(float)
		list_df_annot.append(df_annot)
		if annot_per_geneset == True: # write annot file per annotation per chromosome
			file_out_annot = "{}/{}.{}.{}.annot.gz".format(out_dir, out_prefix, name_annotation, chromosome) # set output filename. ${prefix}.${chr}.annot.gz
			df_annot.to_csv(file_out_annot, sep="\t", index=False, compression="gzip")
		counter += 1
		# if counter == 4: break
	print("CHR={} | Concatenating annotations...".format(chromosome))
	df_annot_combined = pd.concat(list_df_annot, axis='columns') # stack horizontally (there is no joining on indexes, just stacking)
		# *IMPORTANT*: since we did pd.merge(df_bin, df_annot_overlap_bp) with how='left' the know that ALL dfs in list_df_annot have ALL SNPs in df_bim and the order of the SNPs are preserved.
		# ALTERNATIVELY if you don't want 'thin-annot' use this (i.e. adding 'CHR','SNP','CM','BP' columns): df_annot_combined = pd.concat([df_bim]+list_df_annot, axis='columns') # stack horizontally
	
	# print("CHR={} | Calculating standard deviation for annotations...".format(chromosome))
	# df_annot_sd = pd.DataFrame(df_annot_combined.drop(columns=["CHR", "SNP", "CM", "BP"]).std(), columns=["sd"])
	# df_annot_sd.index.name = "annotation"
	# df_annot_sd["n"] = df_annot.shape[1] # number of SNPs in the data frame. This makes it easier to calculate the combined standard deviation across chromosomes later.
	# file_out_annot_combined_sd = "{}/{}.{}.{}.annot_sd".format(args.out_dir, args.out_prefix, "COMBINED_ANNOT", chromosome) 
	# df_annot_sd.to_csv(file_out_annot_combined_sd, sep="\t", index=True)
	### Output file
	### annotation      sd      n
	### antiquewhite3   0.16847050545855485     5
	### blue1   0.1197907423131066      5
	### chocolate       0.0     5

	print("CHR={} | Writing annotations...".format(chromosome))
	if all_genes == True:
		file_out_annot_combined = "{}/all_genes_in_{}.{}.annot.gz".format(out_dir, out_prefix, chromosome)
	else:
		file_out_annot_combined = "{}/{}.{}.{}.annot.gz".format(out_dir, out_prefix, "COMBINED_ANNOT", chromosome)

	df_annot_combined.to_csv(file_out_annot_combined, sep="\t", index=False, compression="gzip")
	
	return None



###################################### MAIN ######################################

out_dir = snakemake.params['out_dir']
out_prefix = snakemake.params['run_prefix']
annot_per_geneset = False
chromosome = snakemake.params['chromosome']
annotations = snakemake.params['annotations']
all_genes = snakemake.params['all_genes']
bimfile = '{}.{}.bim'.format(snakemake.config['LDSC_BFILE_PATH'], chromosome)
dict_of_beds = {}

for name_annot in annotations:
	dict_of_beds[name_annot] = pybedtools.BedTool('{}/{}.{}.bed'.format(out_dir+'/bed', out_prefix, name_annot))

make_annot_file_per_chromosome(chromosome, dict_of_beds, out_dir, out_prefix, annot_per_geneset, bimfile, all_genes)

print("Make annot script is done!")
