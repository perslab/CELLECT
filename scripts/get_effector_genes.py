import pandas as pd
import numpy as np
import os

def human_ens_to_human_symbol(df_unmapped: pd.DataFrame, df_map: pd.DataFrame, drop_unmapped: bool=False, verbose: bool=False) -> None:
    """Map dataframe index human ensemble id to human gene symbol inplace
    Parameters
    ----------
    df_unmapped : DataFrame
        DataFrame with gene names of appropriate format (e.g. ensembl) as index
    
    drop_unmapped : bool, optional (default: False)
        Remove unmapped genes (rows) from df, False: keep original index.
    
    verbose : bool, optional (default: False)
        Print progress report.
    Returns
    -------
        None
    
    TODO
    ----
        * make one mapping-function for all cases
        * modify drop_unmapped to unmapped: {"drop", "keep", "na"}
        * support for custom mapping file
        * handle lower/upper/mixed casee letters in index
    
    """

    assert (len(df_unmapped) > 0), "Empty dataframe."
    
    PREFIX = "ENSG"

    if verbose:
        print("Mapping: human ensembl gene id's --> human gene symbols ...")
    
    # Check that genes are correct format
    mask_peek = np.array([PREFIX in str(idx) for idx in df_unmapped.index.values])

    if not (mask_peek.any()):
        print("Dataframe index contains values that are not ensemble format or not human ensembl id: ", df_unmapped.index.values[mask_peek])
    
#     resource_package = __name__
#     resource_path = 'maps/Homo_sapiens.GRCh38.ens_v90.ensembl2gene_name_version.txt.gz'  # Do not use os.path.join()
#     resource_stream = pkg_resources.resource_stream(resource_package, resource_path)    
#     df_map = pd.read_csv(resource_stream, compression='gzip', delim_whitespace=True)
    # create dictionary for mapping
    map_dict = dict(zip(df_map["ensembl_gene_id"].ravel(), \
                            df_map["gene_name_optimal"].ravel()))

    # map genes in-place,
    # i.e. indexes are replaced directly in df
    df_unmapped.rename(index=map_dict, inplace=True)
    
    if verbose or drop_unmapped:
        # check for unmapped genes
        # note the tilde ~ to get genes NOT mapped
        mask_unmapped = ~df_unmapped.index.isin(df_map["gene_name_optimal"])
        label_unmapped = df_unmapped.index.values[mask_unmapped]
    
        # create report
        n_unmapped = len(label_unmapped)
        
        if verbose:
            n_total = len(df_unmapped)
            pct = n_unmapped / n_total * 100
            print("%.2f pct of genes are unmapped ..." % pct)
        
        if drop_unmapped:
            df_unmapped.drop(index=label_unmapped, inplace=True)
            n_mapped = len(df_unmapped)
            if verbose:
                print("Removed {} unmapped genes ...".format(n_unmapped))
    
    return None



def fnc_get_effector_genes(n_genes_magma, percentile_cutoff_esmu, specificity_id, gwas, es_mu, df_magma, df_map):
    """ take top n magma genes and upper percentile_esmu of non-zero esmu genes per celltype and report intersect
    input args:
    n_genes_magma: number of top magma genes to include (e.g. 100) sorted by p-value in ascending order
    percentile_cutoff_esmu: upper percentile cutoff (e.g. 90) for non-zero ESmu genes to include, sorted by ESmu value in descending order
    specificity_id: expression data id
    es_mu: es_mu table with annotations as columns
    df_magma: <>.resid_correct_all.gsa.genes.out table containing adjusted gene pvals
    df_map: a csv table with column names "ensembl_gene_id", "gene_name_optimal" for mapping genesfrom ensembl to symbol, for passing to subroutine human_ens_to_human_symbol
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
        
        if df_join.shape[0] > 0:
            # append the joined df for the annotation to the others

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

    df_out = df_out.round(decimals={"esmu_percentile":2, "magma_gene_percentile":2, "esmu":2})
    # map genes to symbol 
    df_genes_tmp = pd.DataFrame(data=None, index=df_out["gene_ensembl"])
    human_ens_to_human_symbol(df_unmapped=df_genes_tmp, df_map=df_map, drop_unmapped=False, verbose=True)
    df_out["gene_symbol"] = df_genes_tmp.index
    
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
path_df_map = "data/maps/Homo_sapiens.GRCh38.ens_v90.ensembl2gene_name_version.txt"
print("Finding effector genes for '" + specificity_matrix_name + "' for GWAS '" + gwas_name + "': ")

### Load df_magma
#df_magma = pd.read_csv(base_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".resid_correct_all.gsa.genes.out", sep= '\s+', header = 1)
#NB : FOR TESTING WE USE THE UNCORRECTED PVALS UNTIL CORRECTED ARE AVAILABLE
df_magma = pd.read_csv(magma_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".genes.out", sep= '\s+', header = 0)

### Load esmu df
es_mu = pd.read_csv(specificity_matrix_file, header = 0)

df_map = pd.read_csv(path_df_map, delim_whitespace=True)

df_effector_genes = fnc_get_effector_genes(n_genes_magma=n_genes_magma, 
                                           percentile_cutoff_esmu=percentile_cutoff_esmu, 
                                           specificity_id=specificity_matrix_name, 
                                           gwas=gwas_name, 
                                           es_mu=es_mu, 
                                           df_magma=df_magma,
                                           df_map=df_map)
### Save the results for the given specificity_id and GWAS
outname = specificity_matrix_name + "__" + gwas_name + ".effector_genes.txt"
outdir = cellect_genes_output_dir + "/out"


# dirs should be created by snakemake in fact..
# if not os.path.exists(outdir):              # create the directory manually (to_csv() creates a file if it does not exist, but not a directory)
# 	print("Creating a directory " + outdir)
# 	os.mkdir(outdir) 
# if not os.path.exists(subdir):		    # add the subdirectory
# 	print("Creating a directory " + subdir)
# 	os.mkdir(subdir)	
fullname = os.path.join(outdir, outname)

print("Saving the results...")
df_effector_genes.to_csv(fullname, sep= '\t', index = False)
print("The results are saved as " + fullname)
