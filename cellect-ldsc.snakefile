
# Some overlapping functionality
include: "rules/common_func1.smk"

########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################

# see also the *.smk files

def make_prefix__annotations(prefix, annotations):
	"""
	Makes a list containing the prefix appended to each annotation in the multigeneset file.
	"""
	# This function should possibly be moved into make_cts_file_snake.py - I can't remember
	# why I decided to put it here
	pa_list = []
	for annot in annotations:
		pa_list.append(prefix+'__'+annot)
	return(pa_list)


########################################################################################
################################### VARIABLES ##########################################
########################################################################################

# see also the *.smk files

# Where all the output will be saved
BASE_OUTPUT_DIR = os.path.join(config['BASE_OUTPUT_DIR'], "CELLECT-LDSC")

# More overlapping functionality
include: "rules/common_func2.smk"

SNP_WINDOWS = config['WINDOW_DEFINITION']['WINDOW_LD_BASED']

# Output file prefixes
RUN_PREFIXES = list(SPECIFICITY_INPUT.keys())

wildcard_constraints:
	chromosomes = "\d+",
	BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
	run_prefix = r"|".join(set(SPECIFICITY_INPUT.keys())),
	annotations = r"|".join(set(ANNOTATIONS_DICT)),
	gwas = r"|".join(set(GWAS_SUMSTATS.keys()))


########################################################################################
################################### CONSTANTS ##########################################
########################################################################################

# see also the *.smk files

# These environment variables control how many cores numpy can use
# Setting to 1 allows snakemake to use 1 core per active rule i.e. snakemake core usage = actual core usage
os.environ["MKL_NUM_THREADS"] = str(config['LDSC_CONST']['NUMPY_CORES'])
os.environ["NUMEXPR_NUM_THREADS"] = str(config['LDSC_CONST']['NUMPY_CORES'])
os.environ["OMP_NUM_THREADS"] = str(config['LDSC_CONST']['NUMPY_CORES'])

DATA_DIR = os.path.abspath(config['LDSC_CONST']['DATA_DIR'])
LDSC_DIR = os.path.abspath(config['LDSC_CONST']['LDSC_DIR'])
GENE_COORD_FILE = os.path.abspath(config['GENE_COORD_FILE'])

BFILE_PATH = os.path.join(DATA_DIR,"1000G_EUR_Phase3_plink/1000G.EUR.QC")
PRINT_SNPS_FILE = os.path.join(DATA_DIR,"print_snps.txt")
LD_SCORE_WEIGHTS = os.path.join(DATA_DIR,"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.")
LDSC_BASELINE = os.path.join(DATA_DIR,"baseline_v1.1_thin_annot/baseline.")
SNPSNAP_FILE = os.path.join(DATA_DIR,"ld0.5_collection.tab.gz")
CHROMOSOME_SIZES_FILE = os.path.join(DATA_DIR, "GRCh37-chr-sizes.txt")

SCRIPT_LDSC = os.path.join(LDSC_DIR,'ldsc.py')

CHROMOSOMES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22] 
# OBS: this workflow supports computing LDSC scores for the chromosomes specified in this list. 
# 	   but due to the LDSC software (_N_CHR variable) the LDSC regression have to be run on all chromosomes to work.
# 	   hence the rule 'run_gwas' will fail if not running on all chromosomes.


H2_INTERVAL_ARG_DICT = { # key=mode/out_suffix ; value=cmd_argument
	"qfixed":"--fixed-quantiles",
	### THE BELOW MODES ARE TESTED AND WORKS. They are commented out to simplify the output as many users will not need the files.
	# "q5_exclude_zero":"--exclude0",
	# "q5_with_zero":""
	}


########################################################################################
############################# Pre-check of inputs #######################################
########################################################################################

# see also the *.smk files

if not (config['ANALYSIS_TYPE']['prioritization'] or config['ANALYSIS_TYPE']['conditional'] or config['ANALYSIS_TYPE']['heritability']):
	raise Exception("At least one ANALYSIS_TYPE must be true.")

if config['ANALYSIS_TYPE']['heritability']:
	config_section_string = 'HERITABILITY_INPUT'
	check_conditional_and_heritability_config_sections(config_section_string)
	HERITABILITY_INPUT = build_dict_of_dataset_selected_annotations(config[config_section_string])
	check_conditional_and_heritability_inputs(HERITABILITY_INPUT, ANNOTATIONS_DICT)


if (config['ANALYSIS_TYPE']['heritability_intervals']) and (not config['ANALYSIS_TYPE']['heritability']): 
	raise Exception("Mode 'heritability_intervals' is enabled. This mode requires 'heritability' mode to also be enabled.")


import pandas as pd
# Check GWAS input format for LDSC
def check_gwas_format_for_ldsc(gwas_file):	
	# Column names
	gwas_df = pd.read_csv(gwas_file, sep= '\t')
	if 'N' not in gwas_df.columns:
		raise Exception("Incorrect GWAS input file format: N column is absent: " + gwas_file)
	elif 'SNP' not in gwas_df.columns:
		raise Exception("Incorrect GWAS input file format: SNP column is absent: " + gwas_file)
	elif 'Z' not in gwas_df.columns:
		raise Exception("Incorrect GWAS input file format: Z column is absent: " + gwas_file)

# do it for eqch GWAS in a row
for gwas_name in list(GWAS_SUMSTATS.keys()):
	check_gwas_format_for_ldsc(GWAS_SUMSTATS[gwas_name]['path'])


########################################################################################
################################### Target files ##########################################
########################################################################################

# see also the *.smk files

if config['ANALYSIS_TYPE']['heritability']:
	tmp = "{BASE_OUTPUT_DIR}/results/heritability.csv".format(BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = prefix)
	list_target_files.extend([tmp])
	for prefix in HERITABILITY_INPUT:
		tmp = expand("{BASE_OUTPUT_DIR}/out/h2/{run_prefix}__{gwas}__h2__{annotation}.results",
						run_prefix = prefix,
						BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
						gwas = list(GWAS_SUMSTATS.keys()),
						annotation = HERITABILITY_INPUT[prefix],
						suffix = ["results", "cov", "delete", "part_delete", "log"])
		list_target_files.extend(tmp)
	analysis_types_performed.extend(['heritability'])


if config['ANALYSIS_TYPE']['heritability_intervals']: 
	tmp = "{BASE_OUTPUT_DIR}/results/heritability_intervals.csv".format(BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
	list_target_files.extend([tmp])
	for prefix in HERITABILITY_INPUT:
		tmp = expand('{BASE_OUTPUT_DIR}/out/h2/{run_prefix}__{gwas}__h2_intervals__{annotation}.{mode}.results_intervals',
											run_prefix = prefix,
											BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
											gwas = list(GWAS_SUMSTATS.keys()),
											annotation = HERITABILITY_INPUT[prefix],
											mode=list(H2_INTERVAL_ARG_DICT.keys())),
		list_target_files.extend(tmp)
	analysis_types_performed.extend(['heritability_intervals'])


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
		filter(lambda s: '.csv' not in s, list_target_files) #Not sure if strictly necessary, just to make sure that the .csv files are generated AFTER the analysis
	output:
		expand("{{BASE_OUTPUT_DIR}}/results/{analysis_type}.csv", analysis_type=analysis_types_performed)
	shell:
		"python3 scripts/parse_results.py --base_output_dir {BASE_OUTPUT_DIR}"

rule make_all_genes_background:
	'''
	Makes an all-genes background matrix from each specificity input matrix.
	'''
	input:
		lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path']
	output:
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/all_genes.{run_prefix}.csv"
	conda:
		"envs/cellectpy3.yml"
	params:
		out_dir = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}",
		specificity_matrix_file = lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path'],
		specificity_matrix_name = "{run_prefix}"
	script:
		"scripts/make_all_genes_snake.py"

###################################### CREATE ANNOTATIONS ######################################

if SNP_WINDOWS == True: # Only use SNPs in LD with genes. 

	rule join_snpsnap_bims:
		'''
		Joins SNPsnap file with genes to input BIM file for all chromosomes
		'''
		input:
			"{bfile_path}.CHR_1_22.bim".format(bfile_path = BFILE_PATH),
			"{SNPSNAP_FILE}".format(SNPSNAP_FILE = SNPSNAP_FILE)
		output:
			expand("{{BASE_OUTPUT_DIR}}/precomputation/SNPsnap/SNPs_with_genes.{bfile_prefix}.{chromosome}.txt",
					bfile_prefix = os.path.basename(BFILE_PATH),
					chromosome = CHROMOSOMES) # we don't delete this 'tmp file' because it can be reused for other runs?
		conda:
			"envs/cellectpy3.yml"
		params:
			out_dir = "{BASE_OUTPUT_DIR}/precomputation/SNPsnap",
			chromosomes = CHROMOSOMES,
			snpsnap_file = SNPSNAP_FILE,
			bfile = BFILE_PATH
		script:
			"scripts/join_SNPsnap_and_bim_snake.py"

	rule make_snpsnap_annot:
		'''
		Make the annotation files for input to LDSC from multigeneset files using SNPsnap, LD-based windows
		'''
		input:
			lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path'],
			"{{BASE_OUTPUT_DIR}}/precomputation/SNPsnap/SNPs_with_genes.{bfile_prefix}.{{chromosome}}.txt".format(bfile_prefix = os.path.basename(BFILE_PATH))
		output:
			temp("{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz") # *TEMP FILE*
		conda:
			"envs/cellectpy3.yml"
		params:
			chromosome = "{chromosome}",
			run_prefix = "{run_prefix}",
			precomp_dir = "{BASE_OUTPUT_DIR}/precomputation",
			all_genes = False,
			bfile = BFILE_PATH
		script:
			"scripts/generate_SNPsnap_windows_snake.py"

	rule make_snpsnap_annot_all_genes:
		'''
		Make the annotation files for input to LDSC from multigeneset files using SNPsnap, LD-based windows
		'''
		input:
			"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/all_genes.{run_prefix}.csv",		
			expand("{{BASE_OUTPUT_DIR}}/precomputation/SNPsnap/SNPs_with_genes.{bfile_prefix}.{chromosome}.txt",
					bfile_prefix = os.path.basename(BFILE_PATH),
					chromosome = CHROMOSOMES)
		output:
			"{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.{chromosome}.annot.gz" # NOT TEMP FILE. May be used for h2
		conda:
			"envs/cellectpy3.yml"
		params:
			chromosome = "{chromosome}",
			run_prefix = "{run_prefix}",
			precomp_dir = "{BASE_OUTPUT_DIR}/precomputation",
			all_genes = True,
			bfile = BFILE_PATH
		script:
			"scripts/generate_SNPsnap_windows_snake.py"

else: # Use SNPs in a fixed window size around genes

	rule format_genes:
		'''
		Adds fixed-size windows to either side of all protein-coding genes in the human genome
		'''
		input:
			gene_coords = GENE_COORD_FILE,
			chr_sizes = CHROMOSOME_SIZES_FILE
		output:
			temp(expand("{{BASE_OUTPUT_DIR}}/precomputation/bed/genes_plus_{window}kb.{chromosome}.bed", window = WINDOWSIZE_KB,
																										chromosome = CHROMOSOMES))
		conda:
			"envs/cellectpy3.yml"
		log:
			"{BASE_OUTPUT_DIR}/logs/log.format_genes_snake.txt"
		params:
			windowsize_kb =  WINDOWSIZE_KB,
			bed_out_dir = "{BASE_OUTPUT_DIR}/precomputation/bed"
		script:
			"scripts/format_and_map_snake.py"

	rule find_overlaps:
		'''
		Finds where genes overlap and makes these regions into unique segments in a BED file
		see https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#partition-p-partition
		This is necessary so the max ES score of each overlapping region is only assigned to the region that overlaps
		'''
		input:
			"{BASE_OUTPUT_DIR}/precomputation/bed/genes_plus_{window}kb.{chromosome}.bed"
		output:
			"{BASE_OUTPUT_DIR}/precomputation/bed/overlap_segments_{window}kb.{chromosome}.bed"
		conda:
			"envs/cellectpy3.yml"
		shell:
			"bedops --partition {input} | bedmap --echo --echo-map-id-uniq --delim '\t' - {input} > {output}"
		
	rule make_annot:
		'''
		Make the annotation files used to generate LD scores from the specificity score matrix, HapMap3 SNPs
		and calculated overlapping gene segments
		'''
		input:
			spec_matrix=lambda wildcards: SPECIFICITY_INPUT[wildcards.run_prefix]['path'],
			chrom_bfile="{bfile_path}.{{chromosome}}.bim".format(bfile_path = BFILE_PATH),
			overlap_segments="{{BASE_OUTPUT_DIR}}/precomputation/bed/overlap_segments_{window}kb.{{chromosome}}.bed".format(window = WINDOWSIZE_KB)
		output:
			temp("{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz") # *TEMP FILE*
		log:
			"{BASE_OUTPUT_DIR}/logs/log.make_annot_snake.{run_prefix}.{chromosome}.txt"
		params:
			run_prefix = "{run_prefix}", # better alternative: wildcards.run_prefix?
			chromosome = "{chromosome}",
			out_dir = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}",
			all_genes = False,
			annotations = lambda wildcards: ANNOTATIONS_DICT[wildcards.run_prefix]
		conda:
			"envs/cellectpy3.yml"
		script:
			"scripts/make_annot_from_geneset_all_chr_snake.py"


	rule make_annot_all_genes:
		'''
		Make the annotation files used to generate LD scores from the all-genes matrix, HapMap3 SNPs
		and calculated overlapping gene segments
		'''
		input:
			spec_matrix = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/all_genes.{run_prefix}.csv",
			chrom_bfile = "{bfile_path}.{{chromosome}}.bim".format(bfile_path = BFILE_PATH),
			overlap_segments = "{{BASE_OUTPUT_DIR}}/precomputation/bed/overlap_segments_{window}kb.{{chromosome}}.bed".format(window = WINDOWSIZE_KB)
		output:
			"{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.{chromosome}.annot.gz" # not temp because used in regression
		log:
			"{BASE_OUTPUT_DIR}/logs/log.make_annot_snake.all_genes_in_dataset.{run_prefix}.{chromosome}.txt"
		params:
			run_prefix = "{run_prefix}", # better alternative: wildcards.run_prefix?
			chromosome = "{chromosome}",
			out_dir = "{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset",
			all_genes = True,
			annotations = ["all_genes_in_dataset"]
		conda:
			"envs/cellectpy3.yml"
		script:
			"scripts/make_annot_from_geneset_all_chr_snake.py"


###################################### COMPUTE LDSC SCORES ######################################

rule compute_LD_scores: 
	'''
	Compute the LD scores prior to running LD score regression
	'''
	input:
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz"
	output:
		# ALL these files are tmp files, but it may be an advantage to keep them during pipeline dev, to avoid having to recompting ldscores if something in per_annot changes
		temp("{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.ldscore.gz"),  # *TEMP FILE*
		temp("{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M"), # *TEMP FILE*
		temp("{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M_5_50"), # *TEMP FILE*
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.log" # *TEMP FILE BUT KEEP* 
	wildcard_constraints:
		chromosome="\d+" # chromosome must be only a number, not sure if redundant (also have placed it in this rule arbitrarily)
	log:
		"{BASE_OUTPUT_DIR}/logs/log.compute_LD_scores.{run_prefix}.{chromosome}.txt"
	params:
		chromosome = '{chromosome}',
		run_prefix = '{run_prefix}'
	conda: # Need python 2 for LDSC
		"envs/cellectpy27.yml"
	shell: 
		"{SCRIPT_LDSC} --l2 --bfile {BFILE_PATH}.{params.chromosome} --ld-wind-cm 1 \
		--annot {BASE_OUTPUT_DIR}/precomputation/{params.run_prefix}/{params.run_prefix}.COMBINED_ANNOT.{params.chromosome}.annot.gz \
		--thin-annot --out {BASE_OUTPUT_DIR}/precomputation/{params.run_prefix}/{params.run_prefix}.COMBINED_ANNOT.{params.chromosome} \
		--print-snps {PRINT_SNPS_FILE} &> {log}"

rule compute_LD_scores_all_genes: 
	'''
	Compute the LD scores prior to running LD score regression
	'''
	input:
		"{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.{chromosome}.annot.gz"
	output:
		"{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.{chromosome}.l2.ldscore.gz",
		"{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.{chromosome}.l2.M",
		"{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.{chromosome}.l2.M_5_50",
		"{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.{chromosome}.log"
	wildcard_constraints:
		chromosome="\d+" # chromosome must be only a number, not sure if redundant (also have placed it in this rule arbitrarily)
	log:
		"{BASE_OUTPUT_DIR}/logs/log.compute_LD_scores.all_genes_in_dataset.{run_prefix}.{chromosome}.txt"
	params:
		chromosome = '{chromosome}',
		run_prefix = '{run_prefix}'
	conda: # Need python 2 for LDSC
		"envs/cellectpy27.yml"
	shell: 
		"{SCRIPT_LDSC} --l2 --bfile {BFILE_PATH}.{params.chromosome} --ld-wind-cm 1 \
		--annot {BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{params.run_prefix}.{params.chromosome}.annot.gz \
		--thin-annot --out {BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{params.run_prefix}.{params.chromosome} \
		 --print-snps {PRINT_SNPS_FILE} &> {log}"




for prefix in RUN_PREFIXES:
# Need to use a loop to generate this rule instead of using wildcards because the output depends on the run prefix used.
# The rule generates multiple annotation output files that are matched to a specific run_prefix
# REF https://stackoverflow.com/questions/48993241/varying-known-number-of-outputs-in-snakemake
	rule: # split_LD_scores 
		'''
		Splits the files made during the compute LD scores step by annotation
		'''
		input:
			"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M",
			"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.M_5_50",
			"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.log",
			"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz",
			ldscore="{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/{run_prefix}.COMBINED_ANNOT.{chromosome}.l2.ldscore.gz"
		output:
			expand("{{BASE_OUTPUT_DIR}}/precomputation/{{run_prefix}}/per_annotation/{{run_prefix}}__{annotation}.{{chromosome}}.{suffix}", 
																					annotation=ANNOTATIONS_DICT[prefix], 
																					suffix=["l2.ldscore.gz", "l2.M", "l2.M_5_50", "annot.gz"])
		conda:
			"envs/cellectpy3.yml"
		wildcard_constraints:
			run_prefix = prefix
		params:
			chromosome = '{chromosome}',
			run_prefix = '{run_prefix}',
			out_dir = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}"
		log:
			"{BASE_OUTPUT_DIR}/logs/log.split_ldscores_snake.{run_prefix}_{chromosome}.txt"
		script:
			"scripts/split_ldscores_snake.py"


###################################### PRIORITIZATION + CONDITIONAL ######################################

rule make_cts_file:
	'''
	Makes the cell-type specific file for LDSC cts option
	'''
	input:
		lambda wildcards : expand("{{BASE_OUTPUT_DIR}}/precomputation/{{run_prefix}}/per_annotation/{{run_prefix}}__{annotation}.{chromosome}.l2.ldscore.gz",
																					annotation=ANNOTATIONS_DICT[wildcards.run_prefix],
																					chromosome=CHROMOSOMES)
	output:
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.ldcts.txt"
	conda:
		"envs/cellectpy3.yml"
	params:
		chromosome = CHROMOSOMES,
		prefix__annotations = lambda wildcards: make_prefix__annotations(wildcards.run_prefix, ANNOTATIONS_DICT[wildcards.run_prefix])
	script:
		"scripts/make_cts_file_snake.py"

rule prioritize_annotations:
	'''
	Run LDSC in CTS mode with the provided list of GWAS
	'''
	input:
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.ldcts.txt", # or rule.make_cts_file.output....
		lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['path'],
		expand("{{BASE_OUTPUT_DIR}}/precomputation/control.all_genes_in_dataset/all_genes_in_{{run_prefix}}.{chromosome}.l2.ldscore.gz", 
			chromosome=CHROMOSOMES),
		lambda wildcards: expand("{{BASE_OUTPUT_DIR}}/precomputation/{{run_prefix}}/per_annotation/{{run_prefix}}__{annotation}.{chromosome}.{suffix}",
									annotation=ANNOTATIONS_DICT[wildcards.run_prefix],
									chromosome=CHROMOSOMES,
									suffix=["l2.ldscore.gz", "l2.M", "l2.M_5_50"] # "annot.gz" not needed for CTS mode
									) # files for ALL annotations are listed in the CTS file, so the must be available.
	output:
		"{BASE_OUTPUT_DIR}/out/prioritization/{run_prefix}__{gwas}.cell_type_results.txt"
	log:
		"{BASE_OUTPUT_DIR}/logs/log.prioritize_annotations.{run_prefix}.{gwas}.txt"
	params:
		gwas_path = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['path'], # use wildcards to access dict
		file_out_prefix = '{BASE_OUTPUT_DIR}/out/prioritization/{run_prefix}__{gwas}',
		ldsc_all_genes_ref_ld_chr_name = '{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.'
	conda: # Need python 2 for LDSC
		"envs/cellectpy27.yml"
	shell: 
		"{SCRIPT_LDSC} --h2-cts {params.gwas_path} \
		--ref-ld-chr {LDSC_BASELINE},{params.ldsc_all_genes_ref_ld_chr_name} \
		--w-ld-chr {LD_SCORE_WEIGHTS} \
		--ref-ld-chr-cts {wildcards.BASE_OUTPUT_DIR}/precomputation/{wildcards.run_prefix}.ldcts.txt \
		--out {params.file_out_prefix} &> {log}"

## Conditional
if config['ANALYSIS_TYPE']['conditional']: # needed to ensure CONDITIONAL_INPUT is defined
	rule run_gwas_conditional:
		'''
		Run LDSC in CTS mode conditioned on a single cell-type annotation
		''' 
		input:
			"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.ldcts.txt",
			lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['path'],
			expand("{{BASE_OUTPUT_DIR}}/precomputation/control.all_genes_in_dataset/all_genes_in_{{run_prefix}}.{chromosome}.l2.ldscore.gz", 
																													chromosome=CHROMOSOMES),
			lambda wildcards: expand("{{BASE_OUTPUT_DIR}}/precomputation/{run_prefix}/per_annotation/{run_prefix}__{annotation}.{chromosome}.{suffix}",
									run_prefix=list(CONDITIONAL_INPUT.keys()), # should not be needed in most cases
									annotation=CONDITIONAL_INPUT[wildcards.run_prefix],
									chromosome=CHROMOSOMES,
									suffix=["l2.ldscore.gz", "l2.M", "l2.M_5_50"] # "annot.gz" not needed for CTS mode
									) # files for ALL annotations are listed in the CTS file, so the must be available
		output:
			"{BASE_OUTPUT_DIR}/out/conditional/{run_prefix}__{gwas}__CONDITIONAL__{annotation}.cell_type_results.txt"
		log:
			"{BASE_OUTPUT_DIR}/logs/log.conditional.{run_prefix}.{gwas}.{annotation}.txt"
		params:
			gwas_path = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]['path'],
			file_out_prefix = '{BASE_OUTPUT_DIR}/out/conditional/{run_prefix}__{gwas}__CONDITIONAL__{annotation}',
			ldsc_all_genes_ref_ld_chr_name = "{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.",
			cond_ref_ld_chr_name = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/per_annotation/{run_prefix}__{annotation}."
		conda: # Need python 2 for LDSC
			"envs/cellectpy27.yml"
		shell: 
			"{SCRIPT_LDSC} --h2-cts {params.gwas_path} \
			--ref-ld-chr {LDSC_BASELINE},{params.ldsc_all_genes_ref_ld_chr_name},{params.cond_ref_ld_chr_name} \
			--w-ld-chr {LD_SCORE_WEIGHTS} \
			--ref-ld-chr-cts {wildcards.BASE_OUTPUT_DIR}/precomputation/{wildcards.run_prefix}.ldcts.txt \
			--out {params.file_out_prefix} &> {log}"

########################################################################################
################################### Load rules ##########################################
########################################################################################

if config['ANALYSIS_TYPE']['heritability']: # conditional include statement to speed up building dag. Not sure how effective it is.
	include: "rules/ldsc_h2.smk"



