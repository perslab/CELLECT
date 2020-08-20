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
BASE_OUTPUT_DIR = os.path.join(config['BASE_OUTPUT_DIR'],"CELLECT-MAGMA") # note the overwriting of this variable
CELLECT_GENES_OUTPUT_DIR = os.path.join(config['BASE_OUTPUT_DIR'],"CELLECT-GENES")

# More overlapping functionality
include: "rules/common_func2.smk"

wildcard_constraints:
        CELLECT_GENES_OUTPUT_DIR = CELLECT_GENES_OUTPUT_DIR,
        run_prefix = r"|".join(set(SPECIFICITY_INPUT.keys())),
        annotations = r"|".join(set(ANNOTATIONS_DICT)),
        gwas = r"|".join(set(GWAS_SUMSTATS.keys()))

########################################################################################
############################# Pre-check of inputs ######################################
########################################################################################

if not config['ANALYSIS_TYPE']['effector_genes']:
	raise Exception("config.yml: ANALYSIS TYPE effector_genes must be True.")

import pandas as pd

########################################################################################
################################### PIPELINE ###########################################
########################################################################################



rule all: 
	'''
	Defines the final target files to be generated.
	'''
	input:
		expand("{CELLECT_GENES_OUTPUT_DIR}/out/effector_genes/{run_prefix}__{gwas}.effector_genes.csv", CELLECT_GENES_OUTPUT_DIR=CELLECT_GENES_OUTPUT_DIR, run_prefix=list(SPECIFICITY_INPUT.keys()), gwas=list(GWAS_SUMSTATS.keys()))
        #list_target_files
        
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
                cellect_magma(expand(config['BASE_OUTPUT_DIR']+ "/CELLECT-MAGMA/precomputation/{gwas}/{gwas}.genes.out", gwas=list(GWAS_SUMSTATS.keys()))) #TODO update this to corrected p-values
        output:
                "{CELLECT_GENES_OUTPUT_DIR}/out/effector_genes/{run_prefix}__{gwas}.effector_genes.csv"
        conda:
                "envs/cellectpy3.yml"
        log:
                "{{CELLECT_GENES_OUTPUT_DIR}}/logs/log.get_effector_genes_snake.{run_prefix}.{gwas}.txt"
        params: 
                magma_output_dir = BASE_OUTPUT_DIR,
                cellect_genes_output_dir = CELLECT_GENES_OUTPUT_DIR,
                specificity_matrix_file = lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path'],
                specificity_matrix_name = "{run_prefix}",
                gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id'],
                n_genes_magma = config['N_GENES_MAGMA'],
                percentile_cutoff_esmu = config['PERCENTILE_CUTOFF_ESMU']
        script: 
                "/scripts/get_effector_genes.py"

