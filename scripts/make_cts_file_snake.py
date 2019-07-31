import os

cts_outfile = snakemake.output[0]


prefix__annotations = snakemake.params['prefix__annotations']
chromosome = snakemake.params['chromosome']
file_ldscore = snakemake.input[0] # e.g. <SOMEPATH>/nn_lira_sema.yellow.22.l2.ldscore.gz | output files from make_annot_from_geneset_all_chr.py follow the pattern <PREFIX>.<ANNOTATIONNAME>.<CHR>.l2.ldscore.gz

with open(cts_outfile, "w") as fh_out:
	for pre__annot in prefix__annotations:
		pre__annot_path = "{}/{}.".format(os.path.dirname(file_ldscore), pre__annot) # *OBS*: dot is important
		print(pre__annot + '\t' + pre__annot_path)
		fh_out.write("{}\t{}\n".format(pre__annot, pre__annot_path))

