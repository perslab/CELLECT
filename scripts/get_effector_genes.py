import pandas as pd
import numpy as np
import os

print("effector_genes analysis")

def fnc_get_effector_genes(n_genes_magma, percentile_cutoff_esmu, specificity_id, gwas, es_mu, df_magma):
    """ take top n magma genes and upper percentile_esmu of non-zero esmu genes per celltype and report intersect
    input args:
    n_genes_magma: number of top magma genes to include (e.g. 100) sorted by p-value in ascending order
    percentile_cutoff_esmu: upper percentile cutoff (e.g. 90) for non-zero ESmu genes to include, sorted by ESmu value in descending order
    specificity_id: expression data id
    es_mu: es_mu table with annotations as columns
    df_magma: <>.resid_correct_all.gsa.genes.out table containing adjusted gene pvals
    returns table with columns gwas, specificity_id, annotation, gene_ensembl, gene_symbol, esmu_percentile, magma_gene_percentile, esmu, magma_gene_pval
    """

    # sort df_magma by gene pvals, add percentile column and take top genes
    # TODO: NEEDS UPDATING DEPENDING ON FINAL COLUMN NAME
    df_magma_sorted = df_magma.sort_values(by='P')
    df_magma_sorted["magma_gene_percentile"] = [i/df_magma_sorted.shape[0]*100 for i in reversed(range(df_magma_sorted.shape[0]))]

    df_magma_sorted_head = df_magma_sorted.head(n=n_genes_magma)

    # initialize output df
    df_out = pd.DataFrame(columns=["gwas", "specificity_id", "annotation", "gene_ensembl", "gene_symbol", "esmu_percentile", "magma_gene_percentile", "esmu", "magma_gene_pval"])

    for annotation in es_mu.columns[1:]:
        # remove zero esmu values and sort by esmu
        es_mu_annot_nonzero = es_mu[["gene", annotation]][es_mu[annotation]>0]#.query("`{}`>0".format(annotation))
        es_mu_annot_nonzero.sort_values(by=annotation, ascending=False, inplace=True)
        # add percentile 
        es_mu_annot_nonzero["esmu_percentile"] = [i/es_mu_annot_nonzero.shape[0]*100 for i in reversed(range(es_mu_annot_nonzero.shape[0]))]
        # take the top
        es_mu_annot_nonzero_head = es_mu_annot_nonzero.query("esmu_percentile>{}".format(percentile_cutoff_esmu))
        # join the dataframes on gene
        df_join = pd.merge(es_mu_annot_nonzero_head, df_magma_sorted_head, left_on = "gene", right_on = "GENE", how="inner")
        # append the joined df for the annotation to the others
        print("Got this far with: " + annotation)
        
        df_append = pd.DataFrame({"gwas":[gwas]*df_join.shape[0], 
                                  "specificity_id":[specificity_id]*df_join.shape[0],
                                  "annotation":[annotation]*df_join.shape[0],
                                  "gene_ensembl":df_join["gene"],
                                  "gene_symbol":[None]*df_join.shape[0],
                                  "esmu_percentile":df_join["esmu_percentile"],
                                  "magma_gene_percentile":df_join["magma_gene_percentile"],
                                  "esmu":df_join[annotation],
                                  "magma_gene_pval":df_join["P"]})
        df_out = df_out.append(df_append)

    return df_out


########################################################## MAIN ###############################################################

### Enable logging
snake_log_obj = snakemake.log # class(snakemake.log) = 'snakemake.io.Log
sys.stdout = open(str(snake_log_obj), "w") # could not find better ways than calling str(snake_log_obj) to get the log filepath

### Get the parameters from snakemake
specificity_matrix_name = snakemake.params['specificity_matrix_name']      # ES matrix name
specificity_matrix_file = snakemake.params['specificity_matrix_file']      # ES MU file path
gwas_name = snakemake.params['gwas_name']				   # GWAS name
#base_output_dir = snakemake.params['base_output_dir']			   # base directory for all the outputs
n_genes_magma = snakemake.params['n_genes_magma']
percentile_cutoff_esmu = snakemake.params['percentile_cutoff_esmu']
magma_output_dir = snakemake.params['magma_output_dir']
cellect_genes_output_dir = snakemake.params['cellect_genes_output_dir']

print("Finding effector genes for '" + specificity_matrix_name + "' for GWAS '" + gwas_name + "': ")

### Load df_magma
#df_magma = pd.read_csv(base_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".resid_correct_all.gsa.genes.out", sep= '\s+', header = 1)
#NB : FOR TESTING WE USE THE UNCORRECTED PVALS UNTIL CORRECTED ARE AVAILABLE
df_magma = pd.read_csv(magma_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".genes.out", sep= '\s+', header = 0)

### Load esmu df
es_mu = pd.read_csv(specificity_matrix_file, header = 0)

df_effector_genes = fnc_get_effector_genes(n_genes_magma=n_genes_magma, 
                                           percentile_cutoff_esmu=percentile_cutoff_esmu, 
                                           specificity_id=specificity_matrix_name, 
                                           gwas=gwas_name, 
                                           es_mu=es_mu, 
                                           df_magma=df_magma)
### Save the results for the given specificity_id and GWAS
outname = specificity_matrix_name + "__" + gwas_name + ".effector_genes.csv"
outdir = cellect_genes_output_dir + "/out"
subdir = outdir + "/effector_genes"


# dirs should be created by snakemake in fact..
# if not os.path.exists(outdir):              # create the directory manually (to_csv() creates a file if it does not exist, but not a directory)
# 	print("Creating a directory " + outdir)
# 	os.mkdir(outdir) 
# if not os.path.exists(subdir):		    # add the subdirectory
# 	print("Creating a directory " + subdir)
# 	os.mkdir(subdir)	
fullname = os.path.join(subdir, outname)

print("Saving the results...")
df_effector_genes.to_csv(fullname, index = False)
print("The results are saved as " + fullname)
