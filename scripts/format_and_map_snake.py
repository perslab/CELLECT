import pandas as pd
import numpy as np
import os
import pybedtools

def format_multi_gene_set_file(file_multi_gene_set, out_dir, out_prefix, print_log_files=True):
	""" 
	Reads file_multi_gene_set and formats appropriately based on the argument flags 

	Input
		file_multi_gene_set: file path to a text file (see format specs else where in this file). Supports compressed files  (will be infered from filename, e.g. .gz extension).

	"""
	df_multi_gene_set = pd.read_csv(file_multi_gene_set, sep=None, header=None, engine='python') # sep=None: automatically detect the separator
	df_multi_gene_set.columns = df_multi_gene_set.columns.map(str) # because header=None, the .columns is of type integer (Int64Index(.., dtype='int64')). We need to map array to string before we can rename the columns below
	df_multi_gene_set.columns.values[[0,1,2]] = ["annotation", "gene_input", "annotation_value"]
	if not np.issubdtype(df_multi_gene_set["annotation_value"].dtype, np.number): # REF: https://stackoverflow.com/a/38185759/6639640
		raise Exception("ERROR: your df_multi_gene_set contains non-numeric annotation values. Will not create annotation files.")
	if (df_multi_gene_set["annotation_value"] < 0).any():
		# raise Exception("ERROR: your df_multi_gene_set contains negative annotation values. Will not create annotation files.")
		print("WARNING: your df_multi_gene_set contains negative annotation values. Is this intentional?")
	
	### Check if annotation names are 'valid'
	bool_invalid_annotation_names = df_multi_gene_set["annotation"].str.contains(r"[\s/]|__",regex=True) # character class of whitespace, forward slash and double underscore
	if bool_invalid_annotation_names.any():
		print("file_multi_gene_set header of invalid annotation names:")
		print(df_multi_gene_set[bool_invalid_annotation_names].head(10))
		raise ValueError("The 'annotation' column file_multi_gene_set contains whitespace, forward slash (/) or double underscore (__). Please remove them before running this script.")
		# df_multi_gene_set["annotation"] = df_multi_gene_set["annotation"].replace(r"\s+", "_",regex=True) 
		# ^ any whitespace in annotation_name column in file_multi_gene_set will be converted to underscore ('_').  
		# ^ This is because LDSC .annot files are read as *whitespace delimted* by the ldsc.py program, so annotation_name with whitespace in the name will make the .l2.ldscore.gz header wrong.
		# df_multi_gene_set["annotation"] = df_multi_gene_set["annotation"].replace(r"/", "-",regex=True) 
		# ^ We need to avoid forward slash (/) in the filenames when files are split per annotation (/per_annot dir). 
		# ^ If forward slashes are not replaced, we would get an error when looking for or writing files such as "celltypes.campbell_lvl1.all.campbell_lvl1.a06.NG2/OPC.ges.21.l2.M" when the annotation name is "a06.NG2/OPC.ges"
	### Check for duplicated entries
	bool_duplicated = df_multi_gene_set.duplicated(subset=["annotation", "gene_input"])
	if bool_duplicated.any():
		raise ValueError("file_multi_gene_set contains duplicated genes within an annotation (i.e. non-unique combinations of 'annotation' and 'gene_input' columns). Fix and rerun.")
	print("Read file_multi_gene_set. Header of the parsed/processed file:")
	print(df_multi_gene_set.head(10))
	print("Annotation value summary stats:")
	df_annot_value_sumstats = df_multi_gene_set.groupby("annotation")["annotation_value"].agg(["mean", "std", "max", "min", "count"])
	print(df_annot_value_sumstats)
	file_out_annot_value_sumstatsstats = "{}/log.{}.make_annotation_value_sumstats.txt".format(out_dir, out_prefix)
	if print_log_files and not os.path.exists(file_out_annot_value_sumstatsstats):
		df_annot_value_sumstats.to_csv(file_out_annot_value_sumstatsstats, sep="\t")

	### Ortholog mapping
	map_from_mouse = True
	if df_multi_gene_set.loc[0,'gene_input'].startswith('ENSG'):
		print('\nFirst row contains human gene name, skipping mapping on the assumption all other genes are human also.')
		map_from_mouse = False

	if map_from_mouse == True:
		df_multi_gene_set = map_ensembl_genes_mouse_to_human(df_multi_gene_set, out_dir, out_prefix, print_log_files) # adds "gene" column
	else:
		df_multi_gene_set["gene"] = df_multi_gene_set["gene_input"] # copy


	print("========================== STATS file_multi_gene_set ====================")
	print("Number of gene sets: {}".format(df_multi_gene_set["annotation"].nunique()))
	print("=========================================================================")
	### re-arrange column orders (to make output consistent). 
	# Here we ensure that "annotation", "gene", "annotation_value" columns are always in the same positions.
	# df_multi_gene_set can have any number of columns after these three columns.
	cols = df_multi_gene_set.columns.tolist()
	cols_first = [cols.pop(cols.index(x)) for x in ["annotation", "gene", "annotation_value"]] # OBS: this modifies cols.
	cols = cols_first + cols # extend list
	df_multi_gene_set = df_multi_gene_set[cols] # ALT df.reindex(columns=cols)
	return df_multi_gene_set


def map_ensembl_genes_mouse_to_human(df, out_dir, out_prefix, mouse_to_gene_mapping_path, do_mapping, print_log_files=True):
	"""
	DESCRIPTION
		function to map mouse ensembl genes to human ortholog genes.
	INPUT
	   df: a dataframe with two columns: "annotation" and "gene_input". "gene_input" column should contain Ensembl mouse gene_input IDs to be mapped.
	OUTPUT
	   df: input dataframe with mapped genes. Mouse genes that could not be mapped are removed.
	   file_mapping_summary: a summary file with mapping stats

	REMARKS
		We assume the file_mapping contains only 1-1 mapping (which is true for gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz).
		Otherwise the .map() function might fail
	"""


	print("\nMapping from mouse to human genes")
	### SNIPPET
	# ensembl_gene_id chromosome_name start_position  end_position    mmusculus_homolog_ensembl_gene  mmusculus_homolog_orthology_confidence
	# ENSG00000138593 15      49280673        49338760        ENSMUSG00000035093      1
	# ENSG00000166351 21      14982498        15013906        ENSMUSG00000095294      0
	# ENSG00000168675 18      13217497        13652754        ENSMUSG00000024544      1
	df_mapping = pd.read_csv(mouse_to_gene_mapping_path, delim_whitespace=True, usecols=["ensembl_gene_id", "mmusculus_homolog_ensembl_gene"], index_col=1) # index is MOUSE genes
	### map ortholog genes
	df["gene"] = df["gene_input"].map(df_mapping["ensembl_gene_id"]) # df["gene_input"] is mouse genes. df_mapping["ensembl_gene_id"] is human.
	# ^ .map() returns NaN values for genes not mapped. 
	### make summary of mapping
	file_out_mapping_stats = "{}/log.{}.make_annotation_mapping_stats.txt".format(out_dir, out_prefix)
	if print_log_files and not os.path.exists(file_out_mapping_stats):
		df_summary = df.groupby("annotation")["gene"].agg({'n_genes_input': lambda x: len(x),
															 'n_genes_output': lambda x: len(x)-sum(pd.isnull(x)), 
															'n_genes_not_mapped' : lambda x: sum(pd.isnull(x)),
															'pct_genes_not_mapped': lambda x: "{:.2f}".format(sum(pd.isnull(x))/float(len(x))*100)})
		# ^ we use pd.isnull() instead of np.isnan() because of the issue described here: https://stackoverflow.com/a/36001191/6639640
		df_summary.sort_values(by=['n_genes_output'], inplace=True) 
		df_summary.to_csv(file_out_mapping_stats, sep="\t")
	### final processing
	df.dropna(axis=0, inplace=True) # remove non-mapped genes
	return df

def multi_gene_sets_to_dict_of_beds(df_multi_gene_set, df_gene_coord, windowsize, tmp_bed_dir, out_dir, out_prefix):
	""" 
	INPUT
		df_multi_gene_set: three columns "annotation", "gene" and "annotation_value". Gene is human Ensembl gene names.
	OUTPUT
		dict_of_beds: returns a dict of beds. Keys are annotation names from df_multi_gene_set.
	"""
	print('Making gene set bed files')
	DIR_TMP_PYBEDTOOLS = tmp_bed_dir
	try:
		os.makedirs(DIR_TMP_PYBEDTOOLS, exist_ok=True)
		pybedtools.set_tempdir(DIR_TMP_PYBEDTOOLS) # You'll need write permissions to this directory, and it needs to already exist.
	except Exception as e:
		print("Caught exception: {}".format(e))
	print("Failed setting pybedtools tempdir to {}. Will use standard tempdir /tmp".format(DIR_TMP_PYBEDTOOLS))
	#n_genes_not_in_gene_coord = np.sum(np.isin(df_multi_gene_set["gene"], df_gene_coord["GENE"], invert=True)) # numpy.isin(element, test_elements). Calculates element in test_elements, broadcasting over element only. Returns a boolean array of the same shape as element that is True where an element of element is in test_elements and False otherwise.
	#if n_genes_not_in_gene_coord > 0:
	#    print("*WARNING*: {} genes in the (mapped) input multi gene set is not found in the gene coordinate file. These genes will be discarded".format(n_genes_not_in_gene_coord))
	for name_annotation, df_group in df_multi_gene_set.groupby("annotation"):
		print("Merging input multi gene set with gene coordinates for annotation = {}".format(name_annotation))
		df = pd.merge(df_gene_coord, df_group, left_on="GENE", right_on="gene", how = "inner")
		df['START'] = np.maximum(0, df['START'] - windowsize)
		df['END'] = df['END'] + windowsize
		list_of_lists = [['chr'+(str(chrom).lstrip('chr')), str(start), str(end), str(name), str(score)] for (chrom,start,end,name,score) in np.array(df[['CHR', 'START', 'END', 'GENE', 'annotation_value']])]
		bed_for_annot = pybedtools.BedTool(list_of_lists).sort().merge(c=[4,5], o=["distinct","max"]) 
		bed_for_annot.saveas('{}/{}.{}.bed'.format(out_dir, out_prefix, name_annotation))
	return None


###################################### MAIN ######################################

mouse_to_gene_mapping = snakemake.config['LDSC_MOUSE_GENE_MAPPING']
gene_coords = snakemake.config['LDSC_GENE_COORD_FILE']


windowsize = snakemake.params['windowsize_kb'] * 1000
bed_out_dir = snakemake.params['bed_out_dir']
out_prefix = snakemake.params['run_prefix']

df_multi_gene_set_human = format_multi_gene_set_file(snakemake.input[0], bed_out_dir, out_prefix)
df_gene_coords = pd.read_csv(gene_coords, delim_whitespace = True)
multi_gene_sets_to_dict_of_beds(df_multi_gene_set_human, df_gene_coords, windowsize, bed_out_dir + '/tmp', bed_out_dir, out_prefix)