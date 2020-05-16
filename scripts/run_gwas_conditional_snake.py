import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.tools.tools as sm_tools
import os


	
def fit_multivar_LM(es_mu, df_magma, specificity_id, cond_annotation):
	'''
	Fit linear model between an ESmu matrix 'specificity_id' with one extra ESmu covariate 'cond_annotation' added to each cell type and MAGMA ZSTAT
	'''
	### Extract cell type annotations	
	annotations = es_mu.columns[1:]
	### Delete an annotation that coincides with the conditional annotation
	### to avoid the multicollearity 
	annotations = annotations[annotations != cond_annotation]  
	
	### inner join magma ZSTATs and ES_mu values
	df_regression = pd.merge(es_mu, df_magma, left_on = 'gene', right_on = 'GENE', how = 'inner')
	df_res = pd.DataFrame(columns = ["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"])
	for annotation in annotations:
		y = df_regression.loc[:, df_regression.columns == 'ZSTAT']       # the dependent variable
		X = df_regression.loc[:, df_regression.columns == annotation]    # the 1st independent variable
		# add the 2nd covariate
		X[cond_annotation] = df_regression.loc[:, df_regression.columns == cond_annotation]
		X = sm_tools.add_constant(X.values)      # adding the intercept manually
		ols = sm.OLS(y, X)
		ols_result = ols.fit()
		# FDR correction
		pval = ols_result.pvalues[1]/2                           # get one-sided p-value instead of the two-sided one
		pval = 1 - pval if ols_result.params[1] < 0 else pval    # compute complementary p-value for negative beta's
		# append the result for the given `annotation`
		df_res = df_res.append({"Name": specificity_id + "__" + annotation, "Coefficient": ols_result.params[1], "Coefficient_std_error": ols_result.bse[1], "Coefficient_P_value": pval}, ignore_index = True)
	# return NA when conditioned on itself
	df_res = df_res.append({"Name": specificity_id + "__" + cond_annotation, "Coefficient": np.nan, "Coefficient_std_error": np.nan, "Coefficient_P_value": np.nan}, ignore_index = True)
	df_res = df_res.sort_values(by = ['Coefficient_P_value'])     # sort by original p-value
	return df_res





########################################################## MAIN ###############################################################

### Enable logging
snake_log_obj = snakemake.log # class(snakemake.log) = 'snakemake.io.Log
sys.stdout = open(str(snake_log_obj), "w") # could not find better ways than calling str(snake_log_obj) to get the log filepath

### Get the parameters from snakemake
specificity_matrix_name = snakemake.params['specificity_matrix_name']      # ES matrix name
specificity_matrix_file = snakemake.params['specificity_matrix_file']      # ES MU file path
gwas_name = snakemake.params['gwas_name']				   # GWAS name
base_output_dir = snakemake.params['base_output_dir']			   # base directory for all the outputs
annotations = snakemake.params['annotation']			           # list of cell types for the conditional analysis

### Load MAGMA ZSTATs
df_magma = pd.read_csv(base_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".resid_correct_all.gsa.genes.out", sep= '\s+', header = 1)

### Expression Specificity Metrics
es_mu = pd.read_csv(specificity_matrix_file, header = 0)

for annotation in annotations:
	print("Fitting linear model between MAGMA ZSTATs and ES matrix '" + specificity_matrix_name + "' for GWAS '" + gwas_name + "' conditioned on the cell type annotation'" + annotation  + "': ")

	### Fit the linear model for the given conditional specificity_id and GWAS and one extra ESmu covariate 
	print("Running regressions...")
	df_res = fit_multivar_LM(es_mu, df_magma, specificity_matrix_name, annotation)

	### Save the results
	outname = specificity_matrix_name + "__" + gwas_name + "__CONDITIONAL__" + annotation + ".cell_type_results.txt"
	outdir = base_output_dir + "/out"
	subdir = outdir + "/conditional"
	if not os.path.exists(outdir):              # create the directory manually (to_csv() creates a file if it does not exist, but not a directory)
		print("Creating a directory " + outdir)
		os.mkdir(outdir) 
	if not os.path.exists(subdir):		    # add the subdirectory
		print("Creating a directory " + subdir)
		os.mkdir(subdir)	
	fullname = os.path.join(subdir, outname)

	print("Saving the results...")
	df_res.to_csv(fullname, sep = '\t', index = False, na_rep='NaN')
	print("The results are saved as " + fullname)
