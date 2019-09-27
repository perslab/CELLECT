import pandas as pd
import os
from glob import glob

"""
Description:
    TBD
"""

base_output_dir = snakemake.params['BASE_OUTPUT_DIR']
results_output_dir = snakemake.params['results_out_dir'] 
analysis_types = snakemake.params['analysis_types_performed']


for analysis_type in analysis_types:
    print('Compiling {} result files...'.format(analysis_type))

    if analysis_type == 'prioritization':
        
        prioritization_combined = pd.DataFrame()
        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', analysis_type, '*cell_type_results.txt'))
        result_files = [f for f in glob(result_path)]
    
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
        prioritization_combined.rename(columns = {'Coefficient': 'tau', 'Coefficient_std_error': 'se', 'Coefficient_P_value': 'pvalue'}, inplace = True)
        prioritization_combined.sort_values(['gwas', 'specificity_id'])
        prioritization_combined.to_csv(os.path.join(results_output_dir, 'prioritization.csv'), index = False)


    if analysis_type == 'conditional':
        
        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', analysis_type, '*cell_type_results.txt'))
        result_files = [f for f in glob(result_path)]
        
        conditional_combined = pd.DataFrame()

        for f in result_files:
            
            conditional = pd.read_csv(f, sep = '\t', header = 0)
            
            f = f.split('/')
            metadata = f[-1].replace('.cell_type_results.txt', '') # tabula_muris__EA3_Lee2018__CONDITIONAL__Brain_Non-Myeloid.oligodendrocyte      
            metadata = metadata.split('__') # ['tabula_muris', 'EA3_Lee2018', 'CONDITIONAL', 'Brain_Non-Myeloid.oligodendrocyte']
    
            conditional['gwas'] = metadata[1] # 'EA3_Lee2018'
            conditional['conditional_annotation'] = metadata[-1] # 'Brain_Non-Myeloid.oligodendrocyte'
            
            conditional_combined = pd.concat([conditional_combined, conditional])
        
        
        specificity_id_annotation = conditional_combined['Name'].str.split('__', expand = True) #https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
        conditional_combined['specificity_id'] = specificity_id_annotation[0] #mousebrain
        conditional_combined['annotation'] = specificity_id_annotation[1] #TEGLU32
        conditional_combined = conditional_combined[['gwas', "specificity_id", 'conditional_annotation', 'annotation', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]
        conditional_combined.rename(columns = {'Coefficient': 'tau', 'Coefficient_std_error': 'se', 'Coefficient_P_value': 'pvalue'}, inplace = True)   
        conditional_combined.sort_values(['gwas', 'specificity_id'])
        conditional_combined.to_csv(os.path.join(results_output_dir, 'conditional.csv'), index = False)


    if analysis_type == 'heritability': 

        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', 'h2', '*.results')) #'*results'))
        result_files = [f for f in glob(result_path)]
        
        h2_combined = pd.DataFrame() 

        for f in result_files:
            
            h2 = pd.read_csv(f, sep = '\t', header = 0)
            h2 = h2.tail(1) #.iloc[-1:] didn't exactly as intended, so we try this alternative method

            f = f.split('/')
            metadata = f[-1].replace('.results', '') # tabula_muris__EA3_Lee2018__CONDITIONAL__Brain_Non-Myeloid.oligodendrocyte      
            metadata = metadata.split('__') # ['tabula_muris', 'EA3_Lee2018', 'CONDITIONAL', 'Brain_Non-Myeloid.oligodendrocyte']
            h2['specificity_id'] = metadata[0]
            h2['gwas'] = metadata[1] # 'EA3_Lee2018'
            h2['annotation'] = metadata[-1] # 'Brain_Non-Myeloid.oligodendrocyte'
           
            h2_combined = pd.concat([h2_combined, h2])
           
        h2_combined = h2_combined[['gwas', 'specificity_id', 'annotation', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error', 'Enrichment' , 'Enrichment_std_error', 'Enrichment_p']]
        h2_combined.rename(columns = {'Enrichment': 'h2_enrichment', 'Enrichment_std_error': 'h2_enrichment_se', 'Enrichment_p': 'h2_enrichment_pvalue'}, inplace = True)
        h2_combined.sort_values(['gwas', 'specificity_id'])
        h2_combined.to_csv(os.path.join(results_output_dir, 'heritability.csv'), index = False)


    if analysis_type == 'heritability_intervals':
    
        h2_int_combined = pd.DataFrame()

        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', 'h2', '*results_intervals'))
        result_files = [f for f in glob(result_path )]
        
        for f in result_files:
            
            h2_int = pd.read_csv(f, sep = '\t', header = 0)
            h2_int['q'] = pd.Series(range(h2_int.shape[0]))
            
            f = f.split('/')
            metadata = f[-1].replace('.results_intervals', '') # tabula_muris__EA3_Lee2018__CONDITIONAL__Brain_Non-Myeloid.oligodendrocyte      
            metadata = metadata.split('__') #['tabula_muris', 'BMI_UKBB_Loh2018', 'h2_intervals', 'Liver.hepatocyte.qfixed']
            
            h2_int['specificity_id'] = metadata[0]
            h2_int['gwas'] = metadata[1] # 'BMI_UKBB_Loh2018'
            annotation = metadata[-1] # 'Liver.hepatocyte.qfixed'
            annotation = annotation.split('.')  # ['Liver', 'hepatocyte' , 'qfixed']
            annotation = '.'.join(annotation[:-1]) # 'Liver.hepatocyte'
            h2_int['annotation'] = annotation
            
            h2_int_combined = pd.concat([h2_int_combined, h2_int])
        
            
        h2_int_combined = h2_int_combined[['gwas', 'specificity_id', 'annotation', 'q', 'h2g', 'h2g_se', 'prop_h2g', 'prop_h2g_se', 'enr', 'enr_se', 'enr_pval']]

        h2_int_combined.to_csv(os.path.join(results_output_dir, 'heritability_intervals.csv'), index = False)

    print('{} result files done compiling.'.format(analysis_type))

