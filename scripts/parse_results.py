import pandas as pd
import os
import sys
from glob import glob

import argparse

import logging


"""
Description: parse output files in /out into combined results files in /results.
"""


def main(BASE_OUTPUT_DIR, logger):
    BASE_OUTPUT_DIR = os.path.abspath(BASE_OUTPUT_DIR)
    if not os.path.exists(BASE_OUTPUT_DIR):
        raise IOError("BASE_OUTPUT_DIR={} does not exists.".format(BASE_OUTPUT_DIR))
    RESULTS_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "results")

    result_path_priori = os.path.abspath(os.path.join(BASE_OUTPUT_DIR, 'out', 'prioritization', '*cell_type_results.txt'))
    result_path_cond = os.path.abspath(os.path.join(BASE_OUTPUT_DIR, 'out', 'conditional', '*cell_type_results.txt'))
    result_path_h2 = os.path.abspath(os.path.join(BASE_OUTPUT_DIR, 'out', 'h2', '*.results'))
    result_path_h2_int = os.path.abspath(os.path.join(BASE_OUTPUT_DIR, 'out', 'h2', '*results_intervals'))

    result_files = glob(result_path_priori)
    if result_files:
        logger.info('Compiling prioritization result files from {} files'.format(len(result_files)))
        prioritization_combined = pd.DataFrame()
        for f in result_files:
            prioritization = pd.read_csv(f, sep = '\t', header = 0)
            f = f.split('/')
            metadata = f[-1].replace('.cell_type_results.txt', '') # modules.mousebrain_bmicelltypes__T1D_Bradfield2011        
            metadata = metadata.split('__') # ['modules.mousebrain_bmicelltypes', 'T1D_Bradfield2011']
            prioritization['gwas'] = metadata[-1] # ['T1D_Bradfield2011']
            prioritization_combined = pd.concat([prioritization_combined, prioritization])
        specificity_id_annotation = prioritization_combined['Name'].str.split('__', expand = True) #https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
        prioritization_combined['specificity_id'] = specificity_id_annotation[0] #mousebrain
        prioritization_combined['annotation'] = specificity_id_annotation[1] #TEGLU32
        prioritization_combined = prioritization_combined[['gwas', "specificity_id", 'annotation', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]
        prioritization_combined.rename(columns = {'Coefficient': 'beta', 'Coefficient_std_error': 'beta_se', 'Coefficient_P_value': 'pvalue'}, inplace = True)
        prioritization_combined.sort_values(['gwas', 'specificity_id'])
        prioritization_combined.to_csv(os.path.join(RESULTS_OUTPUT_DIR, 'prioritization.csv'), index = False)

    ### Conditional
    result_files = glob(result_path_cond)
    if result_files:
        logger.info('Compiling conditional result files from {} files'.format(len(result_files)))
        conditional_combined = pd.DataFrame()
        for f in result_files:
            conditional = pd.read_csv(f, sep = '\t', header = 0)                
            f = f.split('/')
            metadata = f[-1].replace('.cell_type_results.txt', '')
            metadata = metadata.split('__')
            conditional['gwas'] = metadata[1]
            conditional['conditional_annotation'] = metadata[-1]
            conditional_combined = pd.concat([conditional_combined, conditional])
        specificity_id_annotation = conditional_combined['Name'].str.split('__', expand = True)
        conditional_combined['specificity_id'] = specificity_id_annotation[0]
        conditional_combined['annotation'] = specificity_id_annotation[1]
        conditional_combined = conditional_combined[['gwas', "specificity_id", 'conditional_annotation', 'annotation', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]
        conditional_combined.rename(columns = {'Coefficient': 'beta', 'Coefficient_std_error': 'beta_se', 'Coefficient_P_value': 'pvalue'}, inplace = True)   
        conditional_combined.sort_values(['gwas', 'specificity_id'])
        conditional_combined.to_csv(os.path.join(RESULTS_OUTPUT_DIR, 'conditional.csv'), index = False)

    result_files = glob(result_path_h2)
    if result_files:
        logger.info('Compiling heritability result files from {} files'.format(len(result_files)))
        h2_combined = pd.DataFrame() 
        for f in result_files:                
            h2 = pd.read_csv(f, sep = '\t', header = 0)
            h2 = h2.tail(1) #.iloc[-1:] didn't exactly as intended, so we try this alternative method
            f = f.split('/')
            metadata = f[-1].replace('.results', '') 
            metadata = metadata.split('__') 
            h2['specificity_id'] = metadata[0]
            h2['gwas'] = metadata[1]
            h2['annotation'] = metadata[-1] 
            h2_combined = pd.concat([h2_combined, h2])
        h2_combined = h2_combined[['gwas', 'specificity_id', 'annotation', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error', 'Enrichment' , 'Enrichment_std_error', 'Enrichment_p']]
        h2_combined.rename(columns = {'Enrichment': 'h2_enrichment', 'Enrichment_std_error': 'h2_enrichment_se', 'Enrichment_p': 'h2_enrichment_pvalue'}, inplace = True)
        h2_combined.sort_values(['gwas', 'specificity_id'])
        h2_combined.to_csv(os.path.join(RESULTS_OUTPUT_DIR, 'heritability.csv'), index = False)

    result_files = glob(result_path_h2_int)
    if result_files:
        logger.info('Compiling heritability_intervals result files from {} files'.format(len(result_files)))
        h2_int_combined = pd.DataFrame()
        for f in result_files:
            h2_int = pd.read_csv(f, sep = '\t', header = 0)
            h2_int['q'] = pd.Series(range(h2_int.shape[0]))
            f = f.split('/')
            metadata = f[-1].replace('.results_intervals', '')
            metadata = metadata.split('__')
            h2_int['specificity_id'] = metadata[0]
            h2_int['gwas'] = metadata[1]
            annotation = metadata[-1]
            annotation = annotation.split('.')
            annotation = '.'.join(annotation[:-1])
            h2_int['annotation'] = annotation
            h2_int_combined = pd.concat([h2_int_combined, h2_int])
        h2_int_combined = h2_int_combined[['gwas', 'specificity_id', 'annotation', 'q', 'h2g', 'h2g_se', 'prop_h2g', 'prop_h2g_se', 'enr', 'enr_se', 'enr_pval']]
        h2_int_combined.to_csv(os.path.join(RESULTS_OUTPUT_DIR, 'heritability_intervals.csv'), index = False)

    logger.info('Result files done compiling.')
        


parser = argparse.ArgumentParser()
parser.add_argument('--base_output_dir', type=str, help="""File path to base CELLECT-LDSC output dir""")
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("parse_results")
main(args.base_output_dir, logger)



### OLD: works with both snakemake params AND commandline args
# try:
#     BASE_OUTPUT_DIR = snakemake.params['BASE_OUTPUT_DIR']
#     # logger = snakemake.logger # --> does not work.
# except:
#     print("parse_results: did not find snakemake.params['BASE_OUTPUT_DIR'] variable. Will use command line argument")
#     if len(sys.argv)<2:
#         print("No command line argument given. Please supply the BASE_OUTPUT_DIR as command line argument.")
#         sys.exit(1)
#     BASE_OUTPUT_DIR = sys.argv[1]
#     print("Got BASE_OUTPUT_DIR from commandline: {}".format(BASE_OUTPUT_DIR))

