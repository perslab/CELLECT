
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
BASE_OUTPUT_DIR = os.path.abspath(config['BASE_OUTPUT_DIR']['MAGMA'])

# More overlapping functionality
include: "rules/common_func2.smk"

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

# Citation info
# The citation info is stored in README.txt, which is saved in magma/bin and duplicated into /data/magma
if not os.path.exists("magma/bin/README.txt"):
	raise Exception("The README file with citation info does not exist in magma/bin.")
f_readme = open("magma/bin/README.txt", "r")
CITATION_INFO = f_readme.read()
f_readme.close()


########################################################################################
############################# Pre-check of inputs ######################################
########################################################################################

if not (config['ANALYSIS_TYPE']['prioritization'] or config['ANALYSIS_TYPE']['conditional']):
        raise Exception("At least one ANALYSIS_TYPE out of 'prioritization' and 'conditional' must be true.")

if (config['ANALYSIS_TYPE']['heritability'] or config['ANALYSIS_TYPE']['heritability_intervals']):
	warnings.warn("'heritability' and 'heritability_intervals' options are available for CELLECT-LDSC only.")

if (config['WINDOW_DEFINITION']['WINDOW_LD_BASED']):
	warnings.warn("WINDOW_LD_BASED is available for CELLECT-LDSC only.")

# see also the *.smk files

########################################################################################
################################### Target files #######################################
########################################################################################

# see the *.smk files

########################################################################################
################################### PIPELINE ###########################################
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


######################## AUTOMATICALLY UNZIP THE GWAS INPUT FILES FOR MAGMA ##############################
        
'''
In this way CELLECT-MAGMA and CELLECT-LDSC can use the same input GWAS files 
'''
rule unzip_gwas:
	input: 
		gwas_file = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['path']
	output:
		temp(expand("{{BASE_OUTPUT_DIR}}/precomputation/{{gwas}}.sumstats"))
	params:
		gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id']
	shell: 
		"gunzip -c {input.gwas_file} > {BASE_OUTPUT_DIR}/precomputation/{params.gwas_name}.sumstats"
		


###################################### GENE ANALYSIS ON SNP P-VALUES ######################################

rule get_uncorrected_pvals:
	'''
	Calculate gene-based P-vals (uncorrected)
	'''
	input: 
		SNPLOC_FILE,
		ANNOT_FILE,
		expand("{{BASE_OUTPUT_DIR}}/precomputation/{gwas}.sumstats", gwas = list(GWAS_SUMSTATS.keys()))
	output: 
		expand("{{BASE_OUTPUT_DIR}}/precomputation/{{gwas}}/{{gwas}}.genes.{ext}", ext = ["raw", "out"])
	params:
		gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id'],
                gwas_file = lambda wildcards: BASE_OUTPUT_DIR + "/precomputation/" + GWAS_SUMSTATS[wildcards.gwas]['id'] + ".sumstats"
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



###################################### PRIORITIZATION + CONDITIONAL ANALYSIS ########################################

rule prioritize_annotations:
	'''
	Fit the linear model between MAGMA ZSTATs and ESmu with the provided list of GWAS
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


# Conditional
if config['ANALYSIS_TYPE']['conditional']: # needed to ensure CONDITIONAL_INPUT is defined
	rule run_gwas_conditional:
		'''
		Run the linear model between MAGMA ZSTATs and ESmu conditioned on a given cell type.
		(One extra ESmu covariate is added to the regression at each step.)
		'''
		input:
			expand("{{BASE_OUTPUT_DIR}}/precomputation/{gwas}/{gwas}.resid_correct_all_ens.gsa.genes.out", gwas = list(GWAS_SUMSTATS.keys())),    # magma ZSTAT files with Ensembl gene IDs
		output:
			expand("{{BASE_OUTPUT_DIR}}/out/conditional/{{run_prefix}}__{{gwas}}__CONDITIONAL__{{annotation}}.cell_type_results.txt")
		log:
			"{BASE_OUTPUT_DIR}/logs/log.conditional.{run_prefix}.{gwas}.{annotation}.txt"
		params:
                	base_output_dir = "{BASE_OUTPUT_DIR}",
                	specificity_matrix_file = lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path'],
                	specificity_matrix_name = "{run_prefix}",
                	gwas_name = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['id'],
			annotation = lambda wildcards: CONDITIONAL_INPUT[wildcards.run_prefix]
		script:
			"scripts/run_gwas_conditional_snake.py"
