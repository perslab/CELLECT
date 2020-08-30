# -*- coding: utf-8 -*-

# Some overlapping functionality
include: "rules/common_func1.smk"

########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################

# see the *.smk files


########################################################################################
################################### VARIABLES ##########################################
########################################################################################

# Where all the output will be saved
BASE_OUTPUT_DIR = os.path.join(config['BASE_OUTPUT_DIR'],"CELLECT-MAGMA") # note the overwriting of this variable, and that we use it to run CELLECT-MAGMA and the below variable for CELLECT-GENES workflow!
CELLECT_GENES_OUTPUT_DIR = os.path.join(config['BASE_OUTPUT_DIR'],"CELLECT-GENES")

# More overlapping functionality
include: "rules/common_func2.smk"

wildcard_constraints:
        BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
        CELLECT_GENES_OUTPUT_DIR = CELLECT_GENES_OUTPUT_DIR,
        run_prefix = r"|".join(set(SPECIFICITY_INPUT.keys())),
        annotations = r"|".join(set(ANNOTATIONS_DICT)),
        gwas = r"|".join(set(GWAS_SUMSTATS.keys()))


########################################################################################
################################### Target files ##########################################
########################################################################################

# see also the *.smk files
# put this here so it is possible to run other cellect workflows without having to set config['ANALYSIS_TYPE']['effector_genes'] to False


list_target_files = []
analysis_types_performed = []     # this is for parsing and compiling the results files

if config['ANALYSIS_TYPE']['effector_genes']:
        tmp = "{CELLECT_GENES_OUTPUT_DIR}/results/effector_genes.csv".format(CELLECT_GENES_OUTPUT_DIR = CELLECT_GENES_OUTPUT_DIR)
        list_target_files.extend([tmp])
        tmp = expand("{CELLECT_GENES_OUTPUT_DIR}/out/{run_prefix}__{gwas}.effector_genes.txt",
                                CELLECT_GENES_OUTPUT_DIR = CELLECT_GENES_OUTPUT_DIR,
                                run_prefix = list(SPECIFICITY_INPUT.keys()),
                                gwas = list(GWAS_SUMSTATS.keys()))
        list_target_files.extend(tmp)
        analysis_types_performed.extend(['effector_genes'])
        
########################################################################################
############################# Pre-check of inputs ######################################
########################################################################################

if not config['ANALYSIS_TYPE']['effector_genes']:
	raise Exception("config.yml: ANALYSIS TYPE effector_genes must be True.")

#if not type(config['N_GENES_MAGMA']) == "int":
#    raise Exception("config.yml: N_GENES_MAGMA must be a positive integer.")
    
if not config['N_GENES_MAGMA']>0:
	raise Exception("config.yml: N_GENES_MAGMA must be a positive integer.")

#if not type(config['PERCENTILE_CUTOFF_ESMU']) == "int":
#    raise Exception("config.yml: PERCENTILE_CUTOFF_ESMU must be an integer between 0 and 100.")

if config['PERCENTILE_CUTOFF_ESMU']<0 or config['PERCENTILE_CUTOFF_ESMU']>100:
	raise Exception("config.yml: PERCENTILE_CUTOFF_ESMU must be an integer between 0 and 100.")

########################################################################################
################################### MODULES ############################################
########################################################################################

import pandas as pd

########################################################################################
################################### PIPELINE ###########################################
########################################################################################

rule all: 
	'''
	Defines the final target files to be generated.
	'''
	input:
		list_target_files #expand("{CELLECT_GENES_OUTPUT_DIR}/out/{run_prefix}__{gwas}.effector_genes.txt",  CELLECT_GENES_OUTPUT_DIR=CELLECT_GENES_OUTPUT_DIR, run_prefix=list(SPECIFICITY_INPUT.keys()), gwas=list(GWAS_SUMSTATS.keys()))
        
     
rule parse_results:
	"""
	Generates {CELLECT_GENES_OUTPUT_DIR}/results/<analysis_type>.csv by parsing ALL output files in {CELLECT_GENES_OUTPUT_DIR}/out/.
	"""
	input:
		filter(lambda s: '.csv' not in s, list_target_files) #Not sure if strictly necessary, just to make sure that the .csv files are generated AFTER the analysis
	output:
		expand("{{CELLECT_GENES_OUTPUT_DIR}}/results/{analysis_type}.csv", analysis_type=analysis_types_performed)
	shell:
		"python3 scripts/parse_results.py --base_output_dir {CELLECT_GENES_OUTPUT_DIR}"


subworkflow cellect_magma:
    snakefile:
        "cellect-magma.snakefile"


rule get_effector_genes:
        '''
        Extract top percentile of ESmu > 0.0 genes 
        Extract top n MAGMA genes
        Take intersection
        '''        
        input: 
                cellect_magma(expand("{BASE_OUTPUT_DIR}/precomputation/{gwas}/{gwas}.resid_correct_all.gsa.genes.pvals.out", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR,gwas=list(GWAS_SUMSTATS.keys()))) 
        output:
                "{CELLECT_GENES_OUTPUT_DIR}/out/{run_prefix}__{gwas}.effector_genes.txt"
        conda:
                "envs/cellectpy3.yml"
        log:
                "{CELLECT_GENES_OUTPUT_DIR}/logs/log.get_effector_genes_snake.{run_prefix}.{gwas}.txt"
        params: 
                magma_output_dir = BASE_OUTPUT_DIR,
                cellect_genes_output_dir = CELLECT_GENES_OUTPUT_DIR,
                specificity_matrix_file = lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path'],
                specificity_matrix_name = "{run_prefix}",
                gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id'],
                n_genes_magma = config['N_GENES_MAGMA'],
                percentile_cutoff_esmu = config['PERCENTILE_CUTOFF_ESMU']
        script: 
                "scripts/get_effector_genes.py"

