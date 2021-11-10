from snakemake.utils import min_version

import sys
import os
import platform
import re
import csv
import gzip
import warnings

min_version("5.27")



_ALLOWED_ID_PATTERN = "^(?!.*__.*)[a-z][a-z0-9_-]+$"


########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################

def check_safe_id(list_of_strings):
        '''
        Returns False if any string in list_of_strings contains the patterns defined in _ILLEGAL_ID_PATTERN.
        '''
        for val in list_of_strings:
                if re.search(_ALLOWED_ID_PATTERN, val):
                        return True
        return False


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


GWAS_SUMSTATS = build_dict_from_id_filepath_key_value_pairs(config['GWAS_SUMSTATS'])
