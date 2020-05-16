import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.tools.tools as sm_tools
import os



	
def fit_LM(specificity_id, es_mu, df_magma):
	### Fit linear model between es and MAGMA ZSTAT

	### Get cell annotations	
	annotations = es_mu.columns[1:]	
	# inner join magma ZSTATs and ES_mu values
	df_regression = pd.merge(es_mu, df_magma, left_on = 'gene', right_on = 'GENE', how = 'inner')
	df_res = pd.DataFrame(columns = ["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"])
	for annotation in annotations:
		y = df_regression.loc[:, df_regression.columns == 'ZSTAT']       # the dependent variable
		X = df_regression.loc[:, df_regression.columns == annotation]
		X = sm_tools.add_constant(X.values)      # adding the intercept manually
		ols = sm.OLS(y, X)
		ols_result = ols.fit()
		# FDR correction
		pval = ols_result.pvalues[1]/2                           # get one-sided p-value instead of the two-sided one
		pval = 1 - pval if ols_result.params[1] < 0 else pval    # compute complementary p-value for negative beta's
		df_res = df_res.append({"Name": specificity_id + "__" + annotation, "Coefficient": ols_result.params[1], "Coefficient_std_error": ols_result.bse[1], "Coefficient_P_value": pval}, ignore_index = True)
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

print("Fitting linear model between MAGMA ZSTATs and ES matrix '" + specificity_matrix_name + "' for GWAS '" + gwas_name + "': ")

### Load MAGMA ZSTATs
df_magma = pd.read_csv(base_output_dir + "/precomputation/" + gwas_name + "/" + gwas_name + ".resid_correct_all.gsa.genes.out", sep= '\s+', header = 1)

### Expression Specificity Metrics
es_mu = pd.read_csv(specificity_matrix_file, header = 0)

### Fit the linear model
print("Running regressions...")
df_res = fit_LM(specificity_matrix_name, es_mu, df_magma)

### Save the results for the given specificity_id and GWAS
outname = specificity_matrix_name + "__" + gwas_name + ".cell_type_results.txt"
outdir = base_output_dir + "/out"
subdir = outdir + "/prioritization"
if not os.path.exists(outdir):              # create the directory manually (to_csv() creates a file if it does not exist, but not a directory)
	print("Creating a directory " + outdir)
	os.mkdir(outdir) 
if not os.path.exists(subdir):		    # add the subdirectory
	print("Creating a directory " + subdir)
	os.mkdir(subdir)	
fullname = os.path.join(subdir, outname)

print("Saving the results...")
df_res.to_csv(fullname, sep = '\t', index = False)
print("The results are saved as " + fullname)
