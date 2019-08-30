
SCRIPT_LDSC_QUANTILE_PERL = os.path.join(LDSC_DIR,'ContinuousAnnotations/quantile_M_fixed_non_zero_quantiles.pl')
SCRIPT_LDSC_H2_RSCRIPT = os.path.join(LDSC_DIR,'ContinuousAnnotations/quantile_h2g.r')

FRQFILE_PREFIX = os.path.join(DATA_DIR, "1000G_Phase3_frq/1000G.EUR.QC.") # needed for h2 estimation



def is_tool(name):
	"""
	Check whether `name` is on PATH and marked as executable.
	Returns boolean value
	"""
	# REF: https://stackoverflow.com/a/34177358/6639640
	from shutil import which
	return which(name) is not None

#########################################################################################
######################################## H2 #############################################
#########################################################################################

if config['ANALYSIS_TYPE']['heritability']: # 'if statement' needed for HERITABILITY_INPUT to be defined. Consider default initialization of HERITABILITY_INPUT
	rule h2: 
		'''
		Estimates h2 for given sets of annotations 
		We add '--print-delete-vals' to enable downstream heritability interval/quantile estimation
		'''
		input: 
			lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]["path"],
			expand("{{BASE_OUTPUT_DIR}}/precomputation/control.all_genes_in_dataset/all_genes_in_{{run_prefix}}.{chromosome}.l2.ldscore.gz", chromosome=CHROMOSOMES),
			lambda wildcards: expand("{{BASE_OUTPUT_DIR}}/precomputation/{{run_prefix}}/per_annotation/{{run_prefix}}__{annotation}.{chromosome}.{suffix}",
				annotation=HERITABILITY_INPUT[wildcards.run_prefix],
				chromosome=CHROMOSOMES,
				suffix=["l2.ldscore.gz", "l2.M", "l2.M_5_50", "annot.gz"] 
				)
		output:
			expand("{{BASE_OUTPUT_DIR}}/out/h2/{{run_prefix}}__{{gwas}}__h2__{{annotation}}.{suffix}", suffix = ["results", "cov", "delete", "part_delete", "log"])
		log:
			"{BASE_OUTPUT_DIR}/logs/log.ldsc_h2.{run_prefix}.{gwas}.{annotation}.txt"
		params:
			gwas_path = lambda wildcards: GWAS_SUMSTATS[wildcards.gwas]["path"],
			file_out_prefix = "{BASE_OUTPUT_DIR}/out/h2/{run_prefix}__{gwas}__h2__{annotation}",
			ldsc_all_genes_ref_ld_chr_name = "{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.",
			h2_ref_ld_chr_name = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/per_annotation/{run_prefix}__{annotation}."
		conda: # Need python 2 for LDSC
			"../envs/cellectpy27.yml"	
		shell:
			"{SCRIPT_LDSC} --h2 {params.gwas_path} \
			--ref-ld-chr {LDSC_BASELINE},{params.ldsc_all_genes_ref_ld_chr_name},{params.h2_ref_ld_chr_name} \
			--frqfile-chr {FRQFILE_PREFIX} \
			--w-ld-chr {LD_SCORE_WEIGHTS} \
			--overlap-annot \
			--thin-annot \
			--print-cov --print-coefficients --print-delete-vals \
			--out {params.file_out_prefix} &> {log}"



#########################################################################################
################################### H2 INTERVAL #########################################
#########################################################################################


# if config['ANALYSIS_TYPE']['heritability_intervals']: 


rule h2_interval_M: 
	'''
	Use quantile_M_fixed_non_zero_quantiles.pl to generate .M files for heritability interval estimation
	'''
	# NOTE: I'm note sure it is necessary to add the 'all genes' to --ref-annot-chr, 
	# But I think it necessary since quantile_h2g.r relies on the .results and .M. files to have the same annotations (in the same order)
	input:
		expand("{{BASE_OUTPUT_DIR}}/precomputation/{{run_prefix}}/per_annotation/{{run_prefix}}__{{annotation}}.{chromosome}.annot.gz", chromosome=CHROMOSOMES),
		expand("{{BASE_OUTPUT_DIR}}/precomputation/control.all_genes_in_dataset/all_genes_in_{{run_prefix}}.{chromosome}.annot.gz", chromosome=CHROMOSOMES)
	output:
		file_M = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/per_annotation/{run_prefix}__{annotation}.{mode}.q_M"
	params:
		ldsc_all_genes_ref_ld_chr_name = '{BASE_OUTPUT_DIR}/precomputation/control.all_genes_in_dataset/all_genes_in_{run_prefix}.',
		file_annot_prefix = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/per_annotation/{run_prefix}__{annotation}.",
		# ^ file_annot_prefix: this is a file prefix. quantile_M_fixed_non_zero_quantiles will add <CHR>.annot.gz to the file
		flag_dependent_argument = lambda wildcards: H2_INTERVAL_ARG_DICT[wildcards.mode]
		# file_out_log = "{}.log".format(os.path.splitext(output['file_M'])[0]) # NOT TESTED
	log:
		"{BASE_OUTPUT_DIR}/logs/log.quantile_M_fixed_non_zero_quantiles.{run_prefix}.{annotation}.{mode}.txt"
	conda:
		"../envs/cellect_perl.yml"
	shell:
		"perl {SCRIPT_LDSC_QUANTILE_PERL} \
		--ref-annot-chr {LDSC_BASELINE},{params.ldsc_all_genes_ref_ld_chr_name},{params.file_annot_prefix} \
		--frqfile-chr {FRQFILE_PREFIX} \
		--annot-header {wildcards.annotation} \
		--nb-quantile 5 \
		--maf 0.05 \
		--thin-annot \
		{params.flag_dependent_argument} \
		--out {output.file_M} &> {log}"

rule h2_interval_rscript: 
	'''
	Runs quantile_h2g.r using h2 .results file and .M file from perl script
	'''
	# Info: quantile_h2g.r
	# Arg1: the file containing the sum of each annotation by quantile of the continuous annotation (e.g. .q5.M file) ['annotfile' inside the script]
	# Arg2: Specify the prefix of the ldsc h2 result outputs (only .results and .part_delete files needed) ['resultfile' inside the script]
	# Arg3: Specify the output filename ['outfile' inside the script]
	input:
		expand("{{BASE_OUTPUT_DIR}}/out/h2/{{run_prefix}}__{{gwas}}__h2__{{annotation}}.{suffix}", suffix = ["results", "part_delete"]),
		file_M = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}/per_annotation/{run_prefix}__{annotation}.{mode}.q_M", # M file output from quantile_M_fixed_non_zero_quantiles.pl
	output:
		file_out_quantile_h2g = '{BASE_OUTPUT_DIR}/out/h2/{run_prefix}__{gwas}__h2_intervals__{annotation}.{mode}.results_intervals'
	params:
		fileout_prefix_ldsc_h2 = "{BASE_OUTPUT_DIR}/out/h2/{run_prefix}__{gwas}__h2__{annotation}" # .<suffixes> are added by Rscript
	# conda:
	# 	"../envs/cellect_R.yml" # You may suddenly experience problems with ldpaths for conda and/or NFS drives.
	# 	# Error from R: "lib/R/etc/ldpaths: No such file or directory"
	# 	# REF: https:/github.com/conda-forge/r-base-feedstock/issues/67
	# 	# SOLUTION --> We drop the R environment and make it a dependency for users.
	log:
		"{BASE_OUTPUT_DIR}/logs/log.quantile_h2g.{run_prefix}.{gwas}.{annotation}.{mode}.txt" # OBS: currently no information in log
	run:
		if not is_tool("Rscript"):
			raise Exception("Could not find Rscript/R on your PATH. Make sure you have installed R when running h2_interval mode")
		shell("Rscript {SCRIPT_LDSC_H2_RSCRIPT} {input.file_M} {params.fileout_prefix_ldsc_h2} {output.file_out_quantile_h2g} &> {log}")




##################################################################################
###################################### WIKI ######################################
##################################################################################


###################################### OUTPUT --h2 ######################################

### Example
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.cov
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.log
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.part_delete
# celltypes.mousebrain__TEGLU4__BMI_UKBB_Loh2018.results


################## OUTPUT quantile_M_fixed_non_zero_quantiles.pl ##################

### WIKI
### The output .M file contains number of SNPs (or annotation size='sum of annotation' for continuos annotation) for each annotation loading in the --ref-annot-chr argument.


### snippet qfixed | <X>.qfixed.M
### Intervals: [0-0.00001],(0-0.2],(0.2-0.4],(0.4-0.6],(0.6-0.8],(0.8-1]
### column 1: annotation size for 'annotation==0' strata
### column 2-6: annotation sizes for annotation strata with fixed equal sizes
# N       4081518 412332  630022  365058  219945  252284
# base    4081518 412332  630022  365058  219945  252284
# Coding_UCSC.bed 36656   9055    15431   11046   6187    6720
# Coding_UCSC.extend.500.bed      165519  42589   69895   45784   27669   28042
# Conserved_LindbladToh.bed       100200  10548   17478   10365   6648    7874
# ...
# WeakEnhancer_Hoffman.bed        69748   11857   18109   11375   6738    7281
# WeakEnhancer_Hoffman.extend.500.bed     297928  49986   76299   47203   27705   30150
# all_genes_in_dataset.mousebrain 2085845 412332  630022  365058  219945  252284
# mousebrain_all.TEINH12.sem_mean 0       40918.9212836192        192911.49975945 178369.311397958        153360.932629917        230021.339404443

### snippet q5_exclude_zero | <X>.q5_exclude_zero.M
### Annotation quintiles are estimated EXCLUDING zero values.
### column 1: annotation_size in lowest quintiles (of the annotation stratified on)
### column 5: annotation_size in highest quintiles (of the annotation stratified on)
# N       376102  376046  376194  377725  373574
# base    376102  376046  376194  377725  373574
# Coding_UCSC.bed 8153    8873    10093   10983   10337
# Coding_UCSC.extend.500.bed      38308   40848   44197   46535   44091
# Conserved_LindbladToh.bed       9821    10410   10155   10900   11627
# ...
# WeakEnhancer_Hoffman.bed        10642   10621   11515   11574   11008
# WeakEnhancer_Hoffman.extend.500.bed     45196   45151   47584   48031   45381
# all_genes_in_dataset.mousebrain 376102  376046  376194  377725  373574
# mousebrain_all.TEINH12.sem_mean 34073.4491382769        95629.6952472526        139792.10504084 206291.635054323        319795.119994766


### snippet q5_with_zero | ....q5_with_zero.M
### Same as q5_exclude_zero, but without excluding zero annotation values when estimating the quintiles


###################################### OUTPUT quantile_h2g.R ######################################

### This output file has one row for each quantile (starting with lowest values) and column summarizing the heritability explained by each quantile, its enrichment and corresponding standard error and P value.
### rows=quantiles ; columns=heritability estimates
# h2g     h2g_se  prop_h2g        prop_h2g_se     enr     enr_se  enr_pval
# 0.0468467691665829      0.00348227231769801     0.283751125243022       0.017433979820968       1.41876013045356        0.0871701758501825      3.7199370574134e-06
# 0.0404121018424434      0.00179110903328923     0.244776311690852       0.00620168879845983     1.22388201557721        0.0310084555740338      1.27746090080597e-11
# 0.0354534348018723      0.00113298619206756     0.214741639556595       0.00132708163840661     1.07370559107702        0.00663539208285362     2.7158960368993e-27
# 0.0292348058573656      0.0013806156620049      0.177075371596941       0.00708277648692335     0.885378842131785       0.0354139617978289      0.00148371877436259
# 0.0131509795779308      0.00262523459741717     0.0796555519125901      0.0162441807658585      0.398276420747926       0.0812206308043398      4.77095371863701e-12


