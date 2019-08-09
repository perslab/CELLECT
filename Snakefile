from snakemake.utils import min_version

import pandas as pd
import numpy as np
import glob
import os

min_version("5.4")

########################################################################################
################################### VARIABLES ##########################################
########################################################################################

configfile: 'config.yml'

BASE_WORKING_DIR = os.path.join(config['BASE_WORKING_DIR'],os.environ['USER'],'cellect-LDSC')
BFILE_PATH = config['LDSC_BFILE_PATH']
PRINT_SNPS_FILE = config['LDSC_PRINT_SNPS_FILE']
LDSC_SCRIPT = config['LDSC_SCRIPT']
MULTIGENESET_DIR = config['LDSC_MULTIGENESET_DIR']
GENE_COORD_FILE = config['LDSC_GENE_COORD_FILE']
GWAS_DIR = config['LDSC_GWAS_DIR']
GWAS_SUMSTATS = config['LDSC_GWAS_SUMSTATS']
LD_SCORE_WEIGHTS = config['LDSC_LD_SCORE_WEIGHTS']
LDSC_BASELINE = config['LDSC_BASELINE']


CHROMOSOMES = config['CHROMOSOMES']
WINDOWSIZE_KB = config['LDSC_WINDOW_SIZE_KB']
RUN_PREFIX = config['RUN_PREFIXES']
 
DATASET = config['DATASET_NAME']
WGCNA = config['WGCNA']
SNP_WINDOWS = config['LDSC_SNPSNAP_WINDOWS']


PRECOMP_DIR = os.path.join(BASE_WORKING_DIR, 'pre-computation')
OUTPUT_DIR = os.path.join(BASE_WORKING_DIR, 'out')


os.environ["MKL_NUM_THREADS"] = str(config['LDSC_NUMPY_MULTITHREADING'])
os.environ["NUMEXPR_NUM_THREADS"] = str(config['LDSC_NUMPY_MULTITHREADING'])
os.environ["OMP_NUM_THREADS"] = str(config['LDSC_NUMPY_MULTITHREADING'])


########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################

def get_annots(run_prefixes, wgcna):
	"""
	Pulls all the annotations from the first multigeneset file.
	"""
	# Should be changed to allow for different cell types in different multigeneset files.
	prefix = run_prefixes[0]
	file_multi_geneset = '{}/multi_geneset.{}.txt'.format(MULTIGENESET_DIR, prefix)
	if wgcna:
		df_multi_gene_set = pd.read_csv(file_multi_geneset, index_col=False)
		annots = np.unique(df_multi_gene_set['module'])
	else:
		df_multi_gene_set = pd.read_csv(file_multi_geneset, sep="\t", index_col=False)
		annots = np.unique(df_multi_gene_set.iloc[:,0])
	return(list(annots))


def make_prefix__annotations(prefix, annotations):
	"""
	Makes a list containing the prefix appended to each annotation in the multigeneset file.
	"""
	pa_list = []
	for annot in annotations:
		pa_list.append(prefix+'__'+annot)
	return(pa_list)

# def get_all_genes_ref_ld_chr_name(dataset, precomp_dir, SNPsnap_windows, windowsize_kb):
# 	""" 
# 	Function to get the ref_ld_chr_name for 'all genes annotation' for ldsc.py --h2/--h2-cts command
# 	"""
# 	# *IMPORTANT*: ldsc_all_genes_ref_ld_chr_name MUST be full file path PLUS trailing "." and prepended ","
# 	window_suffix = ''
# 	if SNPsnap_windows:
# 		SNP_suffix = '_SNP'
# 	else if windowsize_kb != 100:
# 		windowsuffix = '_' + str(windowsize_kb)
# 	ldsc_all_genes_ref_ld_chr_name = ",{precomp_dir}/control.all_genes_in_dataset{win_suffix}/per_annotation/control.all_genes_in_dataset{win_suffix}__all_genes_in_dataset.{dataset}."\
# 										.format(precomp_dir = precomp_dir,
# 											    win_suffix = window_suffix,
# 												dataset = dataset)
# 	return(ldsc_all_genes_ref_ld_chr_name)


########################################################################################
################################### PIPELINE ##########################################
########################################################################################

# Saves the annotations from the first multigeneset as a variable - this should be made into a dictionary so we
# can have a different annotation set for each multigeneset file used as input
ANNOTATIONS = get_annots(RUN_PREFIX, WGCNA)

ANNOTATIONS_ALL_GENES = get_annots(["all_genes_in_dataset"], WGCNA)

rule all: 
	'''
	Defines the final target files to be generated.
	'''
	input:
		expand("{OUTPUT_DIR}/out.ldsc/{run_prefix}__{gwas}.cell_type_results.txt",
			run_prefix = RUN_PREFIX,
			OUTPUT_DIR = OUTPUT_DIR,
			gwas = GWAS_SUMSTATS),
        expand("{PRECOMP_DIR}/control.all_genes_in_dataset/per_annotation/all_genes_in_dataset__{annotation}.{chromosome}.l2.ldscore.gz", 
            PRECOMP_DIR=PRECOMP_DIR,
            annotation = ANNOTATIONS_ALL_GENES,
            chromosome=CHROMOSOMES)



if SNP_WINDOWS == True: # Only use SNPs in LD with genes

	rule join_snpsnap_bims:
		'''
		Joins SNPsnap file with genes to input BIM file for all chromosomes
		'''
		input:
			expand("{bfile_path}.CHR_1_22.bim",
					bfile_path = BFILE_PATH)
		output:
			expand("{{PRECOMP_DIR}}/SNPsnap/SNPs_with_genes.{bfile_prefix}.{chromosome}.txt",
					bfile_prefix = os.path.basename(BFILE_PATH),
					chromosome = CHROMOSOMES)
		conda:
			"envs/cellectpy3.yml"
		params:
			out_dir = "{PRECOMP_DIR}/SNPsnap",
			chromosomes = CHROMOSOMES
		script:
			"scripts/join_SNPsnap_and_bim_snake.py"

	rule make_snpsnap_annot:
		'''
		Make the annotation files for input to LDSC from multigeneset files using SNPsnap, LD-based windows
		'''
		input:
			expand("{MULTIGENESET_DIR}/multi_geneset.{{run_prefix}}.txt",
			 		MULTIGENESET_DIR = MULTIGENESET_DIR),			
			expand("{{PRECOMP_DIR}}/SNPsnap/SNPs_with_genes.{bfile_prefix}.{chromosome}.txt",
					bfile_prefix = os.path.basename(BFILE_PATH),
					chromosome = CHROMOSOMES)
		output:
			"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz"
		conda:
			"envs/cellectpy3.yml"
		params:
			chromosome = "{chromosome}",
			run_prefix = "{run_prefix}",
			precomp_dir = "{PRECOMP_DIR}",
			all_genes = False
		script:
			"scripts/generate_SNPsnap_windows_snake.py"

	rule make_snpsnap_annot_all_genes:
		'''
		Make the annotation files for input to LDSC from multigeneset files using SNPsnap, LD-based windows
		'''
		input:
			expand("{MULTIGENESET_DIR}/multi_geneset.all_genes_in_dataset.txt",
			 		MULTIGENESET_DIR = MULTIGENESET_DIR),			
			expand("{{PRECOMP_DIR}}/SNPsnap/SNPs_with_genes.{bfile_prefix}.{chromosome}.txt",
					bfile_prefix = os.path.basename(BFILE_PATH),
					chromosome = CHROMOSOMES)
		output:
			"{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.annot.gz"
		conda:
			"envs/cellectpy3.yml"
		params:
			chromosome = "{chromosome}",
			run_prefix = "{run_prefix}",
			precomp_dir = "{PRECOMP_DIR}",
			all_genes = True
		script:
			"scripts/generate_SNPsnap_windows_snake.py"

else: # Use SNPs in a fixed window size around genes

	rule format_and_map_genes:
		'''
		Read the multigeneset file, parse and make bed files for each annotation geneset
		'''
		input:
			expand("{MULTIGENESET_DIR}/multi_geneset.{{run_prefix}}.txt",
			 MULTIGENESET_DIR = MULTIGENESET_DIR)
		output:
			expand("{{PRECOMP_DIR}}/{{run_prefix}}/bed/{{run_prefix}}.{annotation}.bed",
					annotation = ANNOTATIONS)
		conda:
			"envs/cellectpy3.yml"
		params:
			wgcna =  WGCNA,
			run_prefix = "{run_prefix}",
			windowsize_kb =  WINDOWSIZE_KB,
			bed_out_dir =  "{PRECOMP_DIR}/{run_prefix}/bed"
		script:
			"scripts/format_and_map_snake.py"

	rule format_and_map_all_genes:
	    '''
	    Works exactly the same way as format_and_map_genes, 
	    but this version was a workaround to overcome
	    the awkward wildcards and to make snakemake 
	    run the same rule twice - on our dataset of interest (fx tabula muris)
	    and on the (control) all_genes_in_dataset
	    '''
	    input:
	        expand("{MULTIGENESET_DIR}/multi_geneset.all_genes_in_dataset.txt",
	                MULTIGENESET_DIR = MULTIGENESET_DIR)
	    output:
	        expand("{{PRECOMP_DIR}}/control.all_genes_in_dataset/bed/all_genes_in_dataset.{annotation}.bed",
	            annotation = ANNOTATIONS_ALL_GENES)  
	    conda:
	        "envs/cellectpy3.yml"
	    params:
	        wgcna =  WGCNA,
	        run_prefix = "all_genes_in_dataset",
	        windowsize_kb =  WINDOWSIZE_KB,
	        bed_out_dir =  "{PRECOMP_DIR}/control.all_genes_in_dataset/bed"
	    script:
	        "scripts/format_and_map_snake.py"

	rule make_annot:
		'''
		Make the annotation files fit for input to LDSC from multigeneset files
		'''
		input:
			expand("{{PRECOMP_DIR}}/{{run_prefix}}/bed/{{run_prefix}}.{annotation}.bed",
					annotation = ANNOTATIONS),
			expand("{bfile_path}.{chromosome}.bim",
					bfile_path = BFILE_PATH,
					chromosome = CHROMOSOMES)
		output:
			"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz"
		conda:
			"envs/cellectpy3.yml"
		params:
			run_prefix = "{run_prefix}",
			chromosome = "{chromosome}",
			out_dir = "{PRECOMP_DIR}/{run_prefix}",
			annotations = ANNOTATIONS
		script:
			"scripts/make_annot_from_geneset_all_chr_snake.py"

	rule make_annot_all_genes:
	    '''
	    Make the annotation files fit for input to to LDSC from multigeneset files
	    '''
	    input:
	        expand("{PRECOMP_DIR}/control.all_genes_in_dataset/bed/all_genes_in_dataset.{annotation}.bed",
	                PRECOMP_DIR = PRECOMP_DIR,
	                annotation = ANNOTATIONS_ALL_GENES),
	        expand("{bfile_prefix}.{chromosome}.bim",
	                bfile_prefix = BFILE_PATH,
	                chromosome = CHROMOSOMES)
	    output:
	        "{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.annot.gz"
	    conda:
	        "envs/cellectpy3.yml"
	    params:
	        run_prefix = "all_genes_in_dataset",
	        chromosome = "{chromosome}",
	        out_dir = PRECOMP_DIR + "/control.all_genes_in_dataset",
	        annotations = ANNOTATIONS_ALL_GENES
	    script:
	        "scripts/make_annot_from_geneset_all_chr_snake.py"


rule compute_LD_scores: 
	'''
	Computing the LD scores prior to running LDSC.
	'''
	input:
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz"
	output:
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.ldscore.gz",
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M",
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M_5_50",
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.log"
	wildcard_constraints:
		chromosome="\d+" # chromosome must be only a number, not sure if redundant (also have placed it in this rule arbitrarily)
	params:
		chromosome = '{chromosome}',
		run_prefix = '{run_prefix}'
	conda: # Need python 2 for LDSC
		"envs/cellectpy27.yml"
	shell: 
		"{LDSC_SCRIPT} --l2 --bfile {BFILE_PATH}.{params.chromosome} --ld-wind-cm 1 --annot {PRECOMP_DIR}/{params.run_prefix}/{params.run_prefix}.COMBINED_ANNOT.{params.chromosome}.annot.gz --thin-annot --out {PRECOMP_DIR}/{params.run_prefix}/{params.run_prefix}.COMBINED_ANNOT.{params.chromosome} --print-snps {PRINT_SNPS_FILE}"

rule compute_LD_scores_all: 
	'''
	Computing the LD scores prior to running LDSC.
	'''
	input:
		"{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.annot.gz"
	output:
		"{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.l2.ldscore.gz",
		"{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.l2.M",
		"{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.l2.M_5_50",
		"{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.log"
	wildcard_constraints:
		chromosome="\d+" # chromosome must be only a number, not sure if redundant (also have placed it in this rule arbitrarily)
	params:
		chromosome = '{chromosome}'
	conda: # Need python 2 for LDSC
		"envs/cellectpy27.yml"
	shell: 
		"{LDSC_SCRIPT} --l2 --bfile {BFILE_PATH}.{params.chromosome} --ld-wind-cm 1 --annot {PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{params.chromosome}.annot.gz --thin-annot --out {PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{params.chromosome} --print-snps {PRINT_SNPS_FILE}"



rule split_LD_scores: 
	'''
	Splits the previously made LD score files by annotation.
	'''
	input:
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.ldscore.gz",
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M",
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M_5_50",
		"{PRECOMP_DIR}/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.log"
	output:
		expand("{{PRECOMP_DIR}}/{{run_prefix}}/per_annotation/{{run_prefix}}__{annotation}.{{chromosome}}.l2.ldscore.gz", annotation=ANNOTATIONS)
	conda:
		"envs/cellectpy3.yml"
	params:
		chromosome = '{chromosome}',
		run_prefix = '{run_prefix}',
		out_dir = "{PRECOMP_DIR}/{run_prefix}"
	script:
		"scripts/split_ldscores_snake.py"

rule split_ldscores_all_genes: 
    '''
    Splits the previously made LD score files by annotation.
    '''
    input:
        "{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.l2.ldscore.gz",
        "{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.l2.M",
        "{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.l2.M_5_50",
        "{PRECOMP_DIR}/control.all_genes_in_dataset/all_genes_in_dataset.COMBINED_ANNOT.{chromosome}.log"
    output:
        expand("{{PRECOMP_DIR}}/control.all_genes_in_dataset/per_annotation/all_genes_in_dataset__{annotation}.{{chromosome}}.l2.ldscore.gz",
                annotation = ANNOTATIONS_ALL_GENES)
    conda:
        "envs/cellectpy3.yml"
    params:
        chromosome = '{chromosome}',
        run_prefix = "all_genes_in_dataset",
        out_dir = "{PRECOMP_DIR}/control.all_genes_in_dataset"
    script:
        "scripts/split_ldscores_snake.py"


rule make_cts_file:
	'''
	Makes the cell-type specific file for LDSC cts flag.
	'''
	input:
		expand("{{PRECOMP_DIR}}/{{run_prefix}}/per_annotation/{{run_prefix}}__{annotation}.{chromosome}.l2.ldscore.gz",
			annotation=ANNOTATIONS,
			chromosome=CHROMOSOMES)
	output:
		"{PRECOMP_DIR}/{run_prefix}.ldcts.txt"
	conda:
		"envs/cellectpy3.yml"
	params:
		chromosome = CHROMOSOMES,
		prefix__annotations = lambda wildcards: make_prefix__annotations(wildcards.run_prefix, ANNOTATIONS)
	script:
		"scripts/make_cts_file_snake.py"

# rule make_cts_file:
#     '''
#     Makes the cell-type specific file for LDSC cts flag.
#     '''
#     input:
#         expand("{PRECOMP_DIR}/{run_prefix}/per_annotation/{run_prefix}__{annotation}.{chromosome}.l2.ldscore.gz",
#             PRECOMP_DIR=PRECOMP_DIR,
#             run_prefix=RUN_PREFIX,
#             annotation=ANNOTATIONS,
#             chromosome=CHROMOSOMES)
#     output:
#         "{PRECOMP_DIR}/{run_prefix}.ldcts.txt"
#     conda:
#         "envs/cellectpy3.yml"
#     params:
#         chromosome = CHROMOSOMES,
#         prefix__annotations = lambda wildcards: make_prefix__annotations(wildcards.run_prefix)
#     script:
#         "scripts/make_cts_file_snake.py"

# rule run_gwas:
# 	'''Run LDSC with the provided list of GWAS'''
# 	input:
# 		expand("{PRECOMP_DIR}/{{run_prefix}}.ldcts.txt", PRECOMP_DIR=PRECOMP_DIR),
# 		expand("{gwas_dir}/{{gwas}}.sumstats.gz", gwas_dir = GWAS_DIR)
# 	output:
# 		"{OUTPUT_DIR}/out.ldsc/{run_prefix}__{gwas}.cell_type_results.txt"
# 	params:
# 		gwas = '{gwas}',
# 		run_prefix = '{run_prefix}',
# 		file_out_prefix = '{OUTPUT_DIR}/out.ldsc/{run_prefix}__{gwas}',
# 		ldsc_all_genes_ref_ld_chr_name = get_all_genes_ref_ld_chr_name(expand("{dataset}", dataset=DATASET)[0],
# 																		expand("{precomp_dir}", precomp_dir=PRECOMP_DIR)[0],
# 																		SNP_WINDOWS)
# 	conda: # Need python 2 for LDSC
# 		"envs/cellectpy27.yml"
# 	shell: 
# 		"{LDSC_SCRIPT} --h2-cts {GWAS_DIR}/{params.gwas}.sumstats.gz --ref-ld-chr {LDSC_BASELINE}{params.ldsc_all_genes_ref_ld_chr_name} --w-ld-chr {LD_SCORE_WEIGHTS} --ref-ld-chr-cts {PRECOMP_DIR}/{params.run_prefix}.ldcts.txt --out {params.file_out_prefix}"


rule run_gwas:
    '''Run LDSC with the provided list of GWAS'''
    input:
        expand("{PRECOMP_DIR}/{{run_prefix}}.ldcts.txt", PRECOMP_DIR=PRECOMP_DIR),
        expand("{gwas_dir}/{{gwas}}.sumstats.gz", gwas_dir = GWAS_DIR),
        expand("{PRECOMP_DIR}/control.all_genes_in_dataset/per_annotation/all_genes_in_dataset__all_genes_in_dataset.{dataset}.{chromosome}.l2.ldscore.gz", 
            PRECOMP_DIR=PRECOMP_DIR,
            dataset = DATASET,
            chromosome=CHROMOSOMES)
    output:
        "{OUTPUT_DIR}/out.ldsc/{run_prefix}__{gwas}.cell_type_results.txt"
    params:
        gwas = '{gwas}',
        run_prefix = '{run_prefix}',
        file_out_prefix = '{OUTPUT_DIR}/out.ldsc/{run_prefix}__{gwas}',
        ldsc_all_genes_ref_ld_chr_name = expand(",{PRECOMP_DIR}/control.all_genes_in_dataset/per_annotation/all_genes_in_dataset__all_genes_in_dataset.{dataset}.", 
                                                dataset = DATASET,
                                                PRECOMP_DIR=PRECOMP_DIR)
    conda: # Need python 2 for LDSC
        "envs/cellectpy27.yml"
    shell: 
        "{LDSC_SCRIPT} --h2-cts {GWAS_DIR}/{params.gwas}.sumstats.gz \
        --ref-ld-chr {LDSC_BASELINE}{params.ldsc_all_genes_ref_ld_chr_name} \
        --w-ld-chr {LD_SCORE_WEIGHTS} \
        --ref-ld-chr-cts {PRECOMP_DIR}/{params.run_prefix}.ldcts.txt \
        --out {params.file_out_prefix}"

