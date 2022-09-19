#########################################################################################
#################################### VARIABLES ##########################################
#########################################################################################

WINDOWSIZE_KB = config['WINDOW_DEFINITION']['WINDOW_SIZE_KB']

SPECIFICITY_INPUT = build_dict_from_id_filepath_key_value_pairs(config['SPECIFICITY_INPUT'])

# Reads the first line of each specificity matrix and saves the annotations
#  as lists where the key is the assigned run prefix
ANNOTATIONS_DICT = get_annots(SPECIFICITY_INPUT)



#########################################################################################
############################## Pre-check of inputs ######################################
#########################################################################################

### Check names/ids
if not check_safe_id(list(SPECIFICITY_INPUT.keys())):
        print(list(SPECIFICITY_INPUT.keys()))
        raise Exception("Illegal charecters in SPECIFICITY_INPUT id's. {}".format(_ALLOWED_ID_MSG))
if not check_safe_id(list(GWAS_SUMSTATS.keys())):
        print(list(GWAS_SUMSTATS.keys()))
        raise Exception("Illegal charecters in GWAS SUMSTATS id's. {}".format(_ALLOWED_ID_MSG))
for key in ANNOTATIONS_DICT:
        for annot in ANNOTATIONS_DICT[key]:
                if not check_safe_id(list(annot)):
                        raise Exception("Illegal charecters in SPECIFICITY_INPUT='{}' with annotation name '{}'. {}".format(key, annot, _ALLOWED_ID_MSG))


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
        