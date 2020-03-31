from snakemake.utils import min_version

import sys
import os
import platform
import re
import csv
import gzip
import warnings

min_version("5.4")



_ILLEGAL_ID_PATTERN = r"\s|__|/"


########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################

def check_safe_id(list_of_strings):
        '''
        Returns False if any string in list_of_strings contains the patterns defined in _ILLEGAL_ID_PATTERN.
        '''
        for val in list_of_strings:
                if re.search(_ILLEGAL_ID_PATTERN, val):
                        return False
        return True


def check_conditional_and_heritability_config_sections(config_section_string):
        if not config_section_string in config:
                raise Exception("Error in config file: parameter {} is required but missing from config file".format(config_section_string))
        try:
                for d in config[config_section_string]: # loop over list of dicts
                        assert isinstance(d["id"], str)
                        assert isinstance(d["annotations"], list)
        except:
                raise Exception("Error in config file: parameter {} is not correctly formatted. Fix the config file and rerun the command".format(config_section_string))


def check_conditional_and_heritability_inputs(dict_dataset_annotations, annotations_dict):
        '''
        dict_dataset_annotations: keys are dataset ids, values are list of 'selected' annotations to perform analysis on (e.g. conditional or h2)
        Checks:
                1. check that dict_dataset_annotations['id'] matches SPECIFICITY_INPUT id's
                2. check that dict_dataset_annotations['annotations'] exists in annotations
        Returns False if does not pass checks, otherwise True
	'''
        for key in dict_dataset_annotations:
                if not key in annotations_dict:
                        raise Exception("[dataset id = {}] used for conditional or heritability analysis but the id was not found in the SPECIFICITY_INPUT. Check your config file".format(key))
                # if not all(dict_dataset_annotations[key] in annotations_dict[key]):
                for annotation in dict_dataset_annotations[key]:
                        if annotation not in annotations_dict[key]:
                                raise Exception("[annotation={}] in [dataset id = {}] in conditional or heritability analysis was not found in the annotations of the SPECIFICITY_INPUT. Check your config file".
format(annotation, key))


def get_annots(specificity_input_dict):
        '''
        Pulls all the annotations from each specificity matrix file and saves them into a dictionary.
        '''
        annots_dict = {}
        for key, dictionary in specificity_input_dict.items():
                if dictionary['path'].endswith('.gz'):
                        fh = gzip.open(dictionary['path'], 'rt') # open in text mode
                else:
                        fh = open(dictionary['path'])
                annotations = next(csv.reader(fh))[1:] # [1:] skip first column because it the the 'gene' column
                annots_dict[key] = annotations # key is dataset name
                fh.close()
        return(annots_dict)


def build_dict_from_id_filepath_key_value_pairs(list_of_dicts):
        '''
        Each dict in the list in MUST contain the keys 'id' and 'path'.
        path will be converted to absolute paths.
        Takes the list of dictionaries and makes it into a new dictionary where the keys are the id values from each dictionary and the values are each dictionary
        e.g. [{"id":"a", "value": 1}, {"id":"b","value":2}] ->
        {"a":{"id":"a", "value": 1}, "b":{"id":"b","value":2}}
        '''
        out_dict = {}
        for d in list_of_dicts:
                d['path'] = os.path.abspath(d['path'])
                out_dict[d['id']] = d
        return(out_dict)


def build_dict_of_dataset_selected_annotations(list_of_dicts):
        '''
        list_of_dicts: list of dicts. Each dict in the list in MUST contain the keys 'id' and 'annotations'.
                list_of_dicts[0]['id']: string
                list_of_dicts[0]['annotations']: list
        returns dict[<id>] = [annotations]
        '''
        dataset_annots_dict = {}
        for d in list_of_dicts:
                dataset_annots_dict[d['id']] = d['annotations']
        return(dataset_annots_dict)


#########################################################################################
#################################### VARIABLES ##########################################
#########################################################################################

#### Load config
## We check if --configfile arg is given to avoid confusing behavior when two config files are loaded.
## snakemake executes the 'configfile: 'config.yml' even if another --configfile is given.
## --configfile will only UPDATE the config dict loaded from 'configfile: 'config.yml'.
## This causes problems if some fields are deleted/missing from the --configfile. Then the config.yml and --configfile will be mixed.
try: # check if config file is already loaded from the --configfile parameter
    config['BASE_OUTPUT_DIR'] # *OBS*: needs to be updated if BASE_OUTPUT_DIR changes name in the config file.
except Exception as e:
    snakemake.logger.info("Loading default config file: config.yml")
    configfile: 'config.yml' # snakemake load config object
else:
        snakemake.logger.info("Loaded config file from --configfile argument") # no Exception raise, so run this


# Where all the output will be saved
BASE_OUTPUT_DIR = os.path.abspath(config['BASE_OUTPUT_DIR'])

WINDOWSIZE_KB = config['WINDOW_DEFINITION']['WINDOW_SIZE_KB']

SPECIFICITY_INPUT = build_dict_from_id_filepath_key_value_pairs(config['SPECIFICITY_INPUT'])
GWAS_SUMSTATS = build_dict_from_id_filepath_key_value_pairs(config['GWAS_SUMSTATS'])

# Reads the first line of each specificity matrix and saves the annotations
#  as lists where the key is the assigned run prefix
ANNOTATIONS_DICT = get_annots(SPECIFICITY_INPUT)


#########################################################################################
#################################### CONSTANTS ##########################################
#########################################################################################

# These environment variables control how many cores numpy can use
# Setting to 1 allows snakemake to use 1 core per active rule i.e. snakemake core usage = actual core usage
os.environ["MKL_NUM_THREADS"] = str(config['MAGMA_CONST']['NUMPY_CORES'])
os.environ["NUMEXPR_NUM_THREADS"] = str(config['MAGMA_CONST']['NUMPY_CORES'])
os.environ["OMP_NUM_THREADS"] = str(config['MAGMA_CONST']['NUMPY_CORES'])


#########################################################################################
############################## Pre-check of inputs ######################################
#########################################################################################

### Check names/ids
if not check_safe_id(list(SPECIFICITY_INPUT.keys())):
        raise Exception("Illegal charecters in SPECIFICITY_INPUT id's. Illegal charecters=[{}]".format(_ILLEGAL_ID_PATTERN))
if not check_safe_id(list(GWAS_SUMSTATS.keys())):
        raise Exception("Illegal charecters in GWAS SUMSTATS id's. Illegal charecters=[{}]".format(_ILLEGAL_ID_PATTERN))
for key in ANNOTATIONS_DICT:
        if not check_safe_id(ANNOTATIONS_DICT[key]):
                raise Exception("Illegal charecters in SPECIFICITY_INPUT={} annotation names. Illegal charecters=[{}]".format(key, _ILLEGAL_ID_PATTERN))


if config['ANALYSIS_TYPE']['conditional']:
        config_section_string = 'CONDITIONAL_INPUT'
        check_conditional_and_heritability_config_sections(config_section_string)
        CONDITIONAL_INPUT = build_dict_of_dataset_selected_annotations(config[config_section_string])
        check_conditional_and_heritability_inputs(CONDITIONAL_INPUT, ANNOTATIONS_DICT)


#########################################################################################
#################################### Target files #######################################
#########################################################################################

list_target_files = []
analysis_types_performed = []     # this is for parsing and compiling the results files

if config['ANALYSIS_TYPE']['prioritization']:
        tmp = "{BASE_OUTPUT_DIR}/results/prioritization.csv".format(BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
        list_target_files.extend([tmp])
        tmp = expand("{BASE_OUTPUT_DIR}/out/prioritization/{run_prefix}__{gwas}.cell_type_results.txt",
                                BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
                                run_prefix = list(SPECIFICITY_INPUT.keys()),
                                gwas = list(GWAS_SUMSTATS.keys()))
        list_target_files.extend(tmp)
        analysis_types_performed.extend(['prioritization'])


if config['ANALYSIS_TYPE']['conditional']:
        tmp = "{BASE_OUTPUT_DIR}/results/conditional.csv".format(BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
        list_target_files.extend([tmp])
        for prefix in CONDITIONAL_INPUT:
                tmp = expand("{BASE_OUTPUT_DIR}/out/conditional/{run_prefix}__{gwas}__CONDITIONAL__{annotation_cond}.cell_type_results.txt",
                                                                        run_prefix = prefix,
                                                                        BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
                                                                        gwas = list(GWAS_SUMSTATS.keys()),
                                                                        annotation_cond = CONDITIONAL_INPUT[prefix])
                list_target_files.extend(tmp)
        analysis_types_performed.extend(['conditional'])

