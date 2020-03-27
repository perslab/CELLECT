from snakemake.utils import min_version

import sys
import os
import platform
import re
import csv
import gzip

min_version("5.4")


_ILLEGAL_ID_PATTERN = r"\s|__|/"


########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################

def check_safe_id(list_of_strings):
        '''
        Returns False if any string in list_of_strings contains the patterns defined in _ILLEGAL_ID_PATTERN.
        '''
        for val in list_of_strings:
                if re.search(_ILLEGAL_ID_PATTERN, val):
                        return False
        return True


def build_dict_from_id_filepath_key_value_pairs(list_of_dicts):
        '''
        Each dict in the list in MUST contain the keys 'id' and 'path'.
        path will be converted to absolute paths.
        Takes the list of dictionaries and makes it into a new dictionary where the keys are the id values from each dictionary and the values are each dictionary
        e.g. [{"id":"a", "value": 1}, {"id":"b","value":2}] ->
        {"a":{"id":"a", "value": 1}, "b":{"id":"b","value":2}}
        '''
        out_dict = {}
        for d in list_of_dicts:
                d['path'] = os.path.abspath(d['path'])
                out_dict[d['id']] = d
        return(out_dict)


def get_annots(specificity_input_dict):
        """
        Pulls all the annotations from each specificity matrix file and saves them into a dictionary.
        """
        annots_dict = {}
        for key, dictionary in specificity_input_dict.items():
                if dictionary['path'].endswith('.gz'):
                        fh = gzip.open(dictionary['path'], 'rt') # open in text mode
                else:
                        fh = open(dictionary['path'])
                annotations = next(csv.reader(fh))[1:] # [1:] skip first column because it the the 'gene' column
                annots_dict[key] = annotations # key is dataset name
                fh.close()
        return(annots_dict)


########################################################################################
################################### VARIABLES ##########################################
########################################################################################

### Load config
# We check if --configfile arg is given to avoid confusing behavior when two config files are loaded.
# snakemake executes the 'configfile: 'config-magma.yml'' even if another --configfile is given.
# --configfile will only UPDATE the config dict loaded from 'configfile: 'config-magma.yml'.
# This causes problems if some fields are deleted/missing from the --configfile. Then the config-magma.yml and --configfile will be mixed.
try: # check if config file is already loaded from the --configfile parameter
    config['BASE_OUTPUT_DIR'] # *OBS*: needs to be updated if BASE_OUTPUT_DIR changes name in the config file.
except Exception as e:
    snakemake.logger.info("Loading default config file: config-magma.yml")
    configfile: 'config-magma.yml' # snakemake load config object
else:
        snakemake.logger.info("Loaded config file from --configfile argument") # no Exception raise, so run this


# Where all CELLECT-MAGMA output will be saved
BASE_OUTPUT_DIR = os.path.abspath(config['BASE_OUTPUT_DIR'])

# Detect OS type and load the corresponding MAGMA binary file
usersystem = platform.system()
magma_version_dir = "magma/bin"
if usersystem == 'Linux':      # Linux
	MAGMA_EXEC = os.path.join(magma_version_dir, "static/magma")
elif usersystem == 'Darwin':   # MacOS
        MAGMA_EXEC = os.path.join(magma_version_dir, "mac/magma")
elif usersystem == 'Windows':  # Win
        raise Exception("Windows OS is not supported at the moment.")
else:                          # the value cannot be determined
	raise Exception("Can not determine the user system/OS name.")

WINDOWSIZE_KB = config['WINDOW_DEFINITION']['WINDOW_SIZE_KB']

SPECIFICITY_INPUT = build_dict_from_id_filepath_key_value_pairs(config['SPECIFICITY_INPUT'])
GWAS_SUMSTATS = build_dict_from_id_filepath_key_value_pairs(config['GWAS_SUMSTATS'])

# Reads the first line of each specificity matrix and saves the annotations
# as lists where the key is the assigned run prefix
ANNOTATIONS_DICT = get_annots(SPECIFICITY_INPUT)


wildcard_constraints:
        BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
        run_prefix = r"|".join(set(SPECIFICITY_INPUT.keys())),
        annotations = r"|".join(set(ANNOTATIONS_DICT)),
        gwas = r"|".join(set(GWAS_SUMSTATS.keys()))


########################################################################################
################################### CONSTANTS ##########################################
########################################################################################

# Data directory
DATA_DIR = os.path.abspath(config['MAGMA_CONST']['DATA_DIR'])

# Data for mapping of human entrez gene IDs to human Ensembl gene IDs
MAPPING_FILE = os.path.join(DATA_DIR, "mapping/gene_id_mapping.hsapiens.ensembl_entrez.txt.gz")
# Gene locations file
GENELOC_FILE = os.path.join(DATA_DIR, "gene_loc.NCBI37.3/NCBI37.3.gene.loc")
# SNP locations file
SNPLOC_FILE = os.path.join(DATA_DIR, "g1000_eur/g1000_eur.bim")
BFILE = os.path.splitext(SNPLOC_FILE)[0]

ANNOT_FILE = os.path.join(BASE_OUTPUT_DIR, "precomputation/NCBI37_1kgp_up" + str(WINDOWSIZE_KB) + "kb_down" + str(WINDOWSIZE_KB) + "kb.genes.annot")
DUMMY_COVAR_FILE_NAME = "magma_dummy_gene_covar_file.NCBI37_3.tab"

# These environment variables control how many cores numpy can use
# Setting to 1 allows snakemake to use 1 core per active rule i.e. snakemake core usage = actual core usage
os.environ["MKL_NUM_THREADS"] = str(config['MAGMA_CONST']['NUMPY_CORES'])
os.environ["NUMEXPR_NUM_THREADS"] = str(config['MAGMA_CONST']['NUMPY_CORES'])
os.environ["OMP_NUM_THREADS"] = str(config['MAGMA_CONST']['NUMPY_CORES'])

# Citation info
# The citation info is stored in README.txt, which is saved in magma/bin and duplicated into /data/magma
if not os.path.exists("magma/bin/README.txt"):
	raise Exception("The README file with citation info does not exist in magma/bin.")
f_readme = open("magma/bin/README.txt", "r")
CITATION_INFO = f_readme.read()
f_readme.close()


########################################################################################
############################# Pre-check of inputs #######################################
########################################################################################

if not (config['ANALYSIS_TYPE']['prioritization'] or config['ANALYSIS_TYPE']['conditional'] or config['ANALYSIS_TYPE']['heritability']):
        raise Exception("At least one ANALYSIS_TYPE must be true.")

### Check names/ids
if not check_safe_id(list(SPECIFICITY_INPUT.keys())):
        raise Exception("Illegal charecters in SPECIFICITY_INPUT id's. Illegal charecters=[{}]".format(_ILLEGAL_ID_PATTERN))
if not check_safe_id(list(GWAS_SUMSTATS.keys())):
        raise Exception("Illegal charecters in GWAS SUMSTATS id's. Illegal charecters=[{}]".format(_ILLEGAL_ID_PATTERN))
for key in ANNOTATIONS_DICT:
        if not check_safe_id(ANNOTATIONS_DICT[key]):
                raise Exception("Illegal charecters in SPECIFICITY_INPUT={} annotation names. Illegal charecters=[{}]".format(key, _ILLEGAL_ID_PATTERN))


if not config['ANALYSIS_TYPE']['prioritization']:
	raise Exception("Currently prioritization is the only available analysis type for MAGMA.")



########################################################################################
################################### Target files ##########################################
########################################################################################

list_target_files = []
analysis_types_performed = []     # this is for parsing and compiling the results files

if config['ANALYSIS_TYPE']['prioritization']:
	tmp = "{BASE_OUTPUT_DIR}/results/prioritization.csv".format(BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
	list_target_files.extend([tmp])
	tmp = expand("{BASE_OUTPUT_DIR}/out/prioritization/{run_prefix}__{gwas}.cell_type_results.txt",
                                BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
                                run_prefix = list(SPECIFICITY_INPUT.keys()),
                                gwas = list(GWAS_SUMSTATS.keys()))
	list_target_files.extend(tmp)
	analysis_types_performed.extend(['prioritization'])



########################################################################################
################################### PIPELINE ##########################################
########################################################################################

rule all:
	'''
	Defines the final target files to be generated.
	'''
	input:
		list_target_files

rule parse_results:
	"""
	Generates {BASE_OUTPUT_DIR}/results/<analysis_type>.csv by parsing ALL output files in {BASE_OUTPUT_DIR}/out/.
	"""
	input:
		filter(lambda s: '.csv' not in s, list_target_files)     # Not sure if strictly necessary, just to make sure that the .csv files are generated AFTER the analysis
	output:
		expand("{{BASE_OUTPUT_DIR}}/results/{analysis_type}.csv", analysis_type = analysis_types_performed)
	shell:
		"python3 scripts/parse_results.py --base_output_dir {BASE_OUTPUT_DIR}"


###################################### CREATE ANNOTATIONS ######################################

rule make_annot:
	'''
	Annotates genes (maps SNPs to genes).
	'''
	input:
		SNPLOC_FILE,
		GENELOC_FILE
	output:
		ANNOT_FILE
	shell:
		"echo \"$(cat magma/bin/README.txt)\"; {MAGMA_EXEC} --annotate window = {WINDOWSIZE_KB},{WINDOWSIZE_KB} \
		--snp-loc {SNPLOC_FILE} \
		--gene-loc {GENELOC_FILE} \
		--out {BASE_OUTPUT_DIR}/precomputation/NCBI37_1kgp_up{WINDOWSIZE_KB}kb_down{WINDOWSIZE_KB}kb"


###################################### CREATE A DUMMY COVAR FILE  ######################################

'''
Create a dummy covar file containing ALL genes in MAGMA universe.
We need to do it because MAGMA only outputs gene scores for genes in the covar file.
'''
rule make_dummy_covar_file:
	output:
		"{BASE_OUTPUT_DIR}/precomputation/" + DUMMY_COVAR_FILE_NAME
	log:
		"{BASE_OUTPUT_DIR}/logs/log.make_dummy_covar_file_snake.txt"
	params:
		geneloc_file = GENELOC_FILE,
		dummy_covar_file = "{BASE_OUTPUT_DIR}/precomputation/" + DUMMY_COVAR_FILE_NAME
	script:
		"scripts/make_dummy_covar_file_snake.py"



###################################### GENE ANALYSIS ON SNP P-VALUES ######################################

rule get_uncorrected_pvals:
	'''
	Calculate gene-based P-vals (uncorrected)
	'''
	input: 
		SNPLOC_FILE,
		ANNOT_FILE,
		lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['path']
	output: 
		expand("{{BASE_OUTPUT_DIR}}/precomputation/{{gwas}}/{{gwas}}.genes.{ext}", ext = ["raw", "out"])
	params:
		gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id'],
                gwas_file = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['path']
	shell:
                "echo \"$(cat magma/bin/README.txt)\"; {MAGMA_EXEC} --bfile {BFILE} \
                --gene-annot {ANNOT_FILE} \
                --pval {params.gwas_file} \
		ncol=N \
		use=SNP,PVAL \
                --out {BASE_OUTPUT_DIR}/precomputation/{params.gwas_name}/{params.gwas_name}"



###################################### GENE-SET ANALYSIS ######################################

rule get_corrected_pvals:
	'''
	Calculate gene-based P-vals (corrected) based on the dummy covar file
	'''
	input:
		dummy_covar_file = "{BASE_OUTPUT_DIR}/precomputation/" + DUMMY_COVAR_FILE_NAME,
		genes_raw = expand("{{BASE_OUTPUT_DIR}}/precomputation/{gwas}/{gwas}.genes.raw", gwas = list(GWAS_SUMSTATS.keys()))
	output:
		expand("{{BASE_OUTPUT_DIR}}/precomputation/{{gwas}}/{{gwas}}.resid_correct_all.{ext}", ext = ["gsa.genes.out", "gsa.out"])
	params:
                gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id']
	shell:
                "echo \"$(cat magma/bin/README.txt)\"; {MAGMA_EXEC} --gene-results {BASE_OUTPUT_DIR}/precomputation/{params.gwas_name}/{params.gwas_name}.genes.raw \
                --gene-covar {input.dummy_covar_file} \
                --model correct=all direction-covar=greater \
		--settings abbreviate=0 gene-info \
                --out {BASE_OUTPUT_DIR}/precomputation/{params.gwas_name}/{params.gwas_name}.resid_correct_all"



################################# MAP ENTREZ HUMAN GENE IDs TO ENSEMBL HUMAN GENE IDs #####################################
 
rule map_human_entrez_to_ens:
	input:
                expand("{{BASE_OUTPUT_DIR}}/precomputation/{gwas}/{gwas}.resid_correct_all.gsa.genes.out", gwas = list(GWAS_SUMSTATS.keys()))    # magma ZSTAT files with Entrez gene IDs		
	output:
		expand("{{BASE_OUTPUT_DIR}}/precomputation/{{gwas}}/{{gwas}}.resid_correct_all_ens.gsa.genes.out")  
	log:
                "{BASE_OUTPUT_DIR}/logs/log.map_human_entrez_to_ens_snake.{gwas}.txt"		
	params:
                base_output_dir = "{BASE_OUTPUT_DIR}",
                gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id'],
		mapping_file = MAPPING_FILE
	script:
		"scripts/map_human_entrez_to_ens.py"



###################################### PRIORITIZATION ########################################

rule prioritize_annotations:
	'''
	Fit linear model between MAGMA ZSTATs and ES with the provided list of GWAS
	'''
	input:
		expand("{{BASE_OUTPUT_DIR}}/precomputation/{gwas}/{gwas}.resid_correct_all_ens.gsa.genes.out", gwas = list(GWAS_SUMSTATS.keys())),    # magma ZSTAT files with Ensembl gene IDs
		lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path']      # es mu files
	output:
		expand("{{BASE_OUTPUT_DIR}}/out/prioritization/{{run_prefix}}__{{gwas}}.cell_type_results.txt")
	log:
		"{BASE_OUTPUT_DIR}/logs/log.prioritize_annotations_snake.{run_prefix}.{gwas}.txt"
	params:
                base_output_dir = "{BASE_OUTPUT_DIR}",
                specificity_matrix_file = lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path'],
                specificity_matrix_name = "{run_prefix}",
                gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id']
	script:
		"scripts/prioritize_annotations_snake.py"


