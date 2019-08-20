# CELLECT

**CELL**-type **E**xpression-specific integration for **C**omplex **T**raits (**CELLECT**) is a computational toolkit for  identifing likely etiologic cell-types underlying complex traits. CELLECT leverages existing genetic prioritization models to integrate single-cell transcriptomic and human genetic data when identifing likely etiologic cell-types. Currently we have only implemented CELLECT-LDSC that uses [LDSC](https://github.com/bulik/ldsc) for genetic prioritization.

![fig-integration](https://user-images.githubusercontent.com/5487016/62281981-0cb33d00-b44f-11e9-8c0b-24aaa2b7d286.png)


## Update log 

**August 1st 2019: v0.1**. Beta version implementing CELLECT-LDSC. We expect to release CELLECT-MAGMA, CELLECT-RolyPoly and CELLECT-DEPICT in the near future.

**August 12th 2019: v0.1.1**. CELLECT now takes a matrix as input and CELLECT-LDSC does not require an all_genes background as this is generated from the matrix.

## How does CELLECT work?

CELLECT quantifies the association between common polygenetic GWAS signal (heritability) and cell-type expression specificity (ES) of genes using established genetic prioritization models such as LDSC (Hilary Kiyo Finucane et al., 2015), RolyPoly (Calderon et al., 2017), DEPICT (Pers et al., 2015) or MAGMA covariate analysis (Skene et al., 2018).

CELLECT takes as input GWAS data and cell-type expression specificity estimates. In order to compute robust estimates of ES, we developed the computational method called **CELLEX** (**CELL**-type **EX**pression-specificity). CELLEX is built on the observation that different ES metrics provide complementary cell-type expression specific profiles. Our method incorporates a ‘wisdom of the crowd’ approach by integrating multiple ES metrics to obtain improved robustness and a more expressive ES measure that captures multiple aspects of expression specificity.  CELLEX can be found [here](https://github.com/perslab/CELLEX).

![fig-CELLECT-conceptual-h2](https://user-images.githubusercontent.com/5487016/62367093-e3ff7600-b528-11e9-8879-8f69005fbea5.png)

Schematic illustration of CELLECT and CELLEX. The bottom layer shows a disease with multiple genetic components (G1-G4). CELLECT integrates disease heritability estimates with cell-type expression specificity to identify the etiologic cell-types (T1 and T4) underlying the genetic components (G1 and G4). CELLEX estimates expression specificity from single-cell transcriptomic atlases.


## Installation

1. **Clone CELLECT repository**
    ```
    git clone --recurse-submodules https://github.com/perslab/CELLECT.git
    ```
The `--recurse-submodules` is needed to clone the [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) 'ldsc' ([pascaltimshel/ldsc](https://github.com/pascaltimshel/ldsc)), which is a modfied version of the original ldsc repository.

2. **Install Snakemake via conda**

   CELLECT uses the workflow management software [**Snakemake**](https://snakemake.readthedocs.io/en/stable/). To make things easier for you, CELLECT snakemake workflow utilises **conda environments** to avoid any issues with software dependencies and versioning. CELLECT snakemake workflow will automatically install all necessary dependencies. All you need to do is to [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) (if conda is not already present on your system) and then install snakemake:
   
    ```bash
    conda install -c bioconda snakemake
    ```
3. **Download CELLECT-data**
    
    There are several data files (roughly 6GB after unpacking) which are required to use CELLECT they can be downloaded as follows:
    
    ```bash
    wget www.WILLADDLINKLATER.com/cellect-data.gz
    tar -xvf cellect-data.gz
    ```

## Getting started with CELLECT-LDSC


1. **Modify the `config.yml` file**. The config file is divided up into two categories (plus model-specific sub categories):

* **RUN-SPECIFIC**: These variables can change from run-to-run and affect how CELLECT processes its input data.
* **CONSTANTS**: These variables include things like paths to data and scripts that generally do not change between runs. Thus you typically only need to update them once for each system.


2. When the config file contains the above we can run CELLECT-LDSC by navigating to the cloned CELLECT directory and entering the following:

```bash
snakemake --use-conda -j 4
```
This will run the workflow using 4 cores (`-j 4`). If you wish to use to use all available cores pass just the `-j` flag.


### CELLECT-LDSC Example: 

```
COMING SOONG

cd CELLECT
wget www.ssgac.com/some_sumstats.gz example/XXX.gwassumstats.gz
(munge?)
snakemake --use-conda --configfile -j example/cellect-ldsc.yml
```

```
TODO: add data/example/{...} containing example CELLEX output
TODO: add info for downloading 1 GWAS sum stat. 
TODO: add example/cellect-ldsc.yml config file with preconfigured variables
```

## Documentation

See the below sections for relevant documentation and information about parameters, input and output formats.

### `SPECIFICITY_INPUT`: Expression specificity

This is one or several matrices containing genes in the index/first column and annotations in the subsequent columns. The gene names must be in **Ensembl human** format and the specificity values in should be between **0 and 1** e.g.

| **gene** 			  | **Bladder.bladder_cell**  | ... | **Trachea.mesenchymal_cell** |
|---------------------|---------------------------|-----|------------------------------|
| **ENSG00000081791** | 0.43                      | ... | 0.11                         |
| ...                 | ...                       | ... | ...                          |
| **ENSG00000101440** | 0.21                      | ... | 0.89                         |


In the `config.yml` please provide specificity input as both the name you would like your output files to have and a path to the matrix.

During the LD score regression, the list of genes in the index column (with each gene assigned a value of 1) will be used  to create an 'all genes' multigeneset to use as an indepedent variable.

### `GWAS_SUMSTATS`: GWAS summary statistics

This must be one of several already munged (using the `munge-sumstats.py` script found in LD score regression) summary statistics for a given trait.


In the `config.yml` please provide GWAS summary stats as both the name would like your output files to have and a path to the summary statistic file.

### `BASE_OUTPUT_DIR`: Output directory

This is a path to a directory where you would like your output to be saved. Ideally use a path to a solid state drive to speed up computation. **2-3 GBs of space are usually needed** for each Specificity Input but additional GWAS summary stats will not take up much more storage.



