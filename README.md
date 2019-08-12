# CELLECT

**CELL**-type **E**xpression-specific integration for **C**omplex **T**raits (**CELLECT**) is a computational toolkit for  identifing likely etiologic cell-types underlying complex traits. CELLECT leverages existing genetic prioritization models to integrate single-cell transcriptomic and human genetic data when identifing likely etiologic cell-types. 

![fig-integration](https://user-images.githubusercontent.com/5487016/62281981-0cb33d00-b44f-11e9-8c0b-24aaa2b7d286.png)


### Update log 

**August 1st 2019: v0.1**. Beta version. Currently CELLECT only support S-LDSC software for genetic identification of etiologic cell-types. We expect to release a pipeline for MAGMA, RolyPoly and DEPICT in the near future. 

## How does it work?

CELLECT quantifies the association between common polygenetic GWAS signal (heritability) and cell-type expression specificity (ES) of genes using established genetic prioritization models such as S-LDSC (Hilary Kiyo Finucane et al., 2015), RolyPoly (Calderon et al., 2017), DEPICT (Pers et al., 2015) or MAGMA covariate analysis (Skene et al., 2018).

CELLECT takes as input GWAS data and cell-type expression specificity estimates. In order to compute robust estimates of ES, we developed the computational method called **CELLEX** (**CELL**-type **EX**pression-specificity). CELLEX is built on the observation that different ES metrics provide complementary cell-type expression specific profiles. Our method incorporates a ‘wisdom of the crowd’ approach by integrating multiple ES metrics to obtain improved robustness and a more expressive ES measure that captures multiple aspects of expression specificity.  CELLEX can be found [here](https://github.com/perslab/CELLEX).

![fig-CELLECT-conceptual-h2](https://user-images.githubusercontent.com/5487016/62367093-e3ff7600-b528-11e9-8879-8f69005fbea5.png)

Schematic illustration of CELLECT and CELLEX. The bottom layer shows a disease with multiple genetic components (G1-G4). CELLECT integrates disease heritability estimates with cell-type expression specificity to identify the etiologic cell-types (T1 and T4) underlying the genetic components (G1 and G4). CELLEX estimates expression specificity from single-cell transcriptomic atlases.


## Installation

After cloning this repository there are a few other things that need to be set-up before CELLECT can be run.


Our forked version of this repository [(pascaltimshel/ldsc)](https://github.com/pascaltimshel/ldsc) **must also be cloned** for CELLECT-LDSC to work.

The models in CELLECT are built in **Snakemake** and the pipelines utilise **conda environments**. Therefore the easiest way to get started would be to [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) (if conda is not already present on your system) and then either within base conda or a conda environment:
```bash
conda install -c bioconda snakemake
```

There are several data files (roughly 6GB after unpacking) which are required to use CELLECT they can be downloaded as follows:
```bash
wget www.linktocellect.com/cellect-data.gz
tar -xvf cellect-data.gz
```

Finally, within your cloned version of this repository, modify the `config.yml` file so that it is specific to the system you are working on. The config file is divided up into two categories (plus model-specific sub categories):

* **RUN-SPECIFIC**: These variables can change from run-to-run and affect how CELLECT processes its input data.
* **CONSTANTS**: These will includes things like paths to data and scripts that we do not expect to change between runs.

For now, only modify the CONSTANTS part to point to the cloned LDSC script and the data directory we just downloaded.

## Usage of CELLECT-LDSC

To run CELLECT-LDSC we must specify in the `config.yml` several things:

### Specificity input

This is one or several matrices containing genes in the index/first column and annotations in the subsequent columns. The gene names must be in **Ensembl human** format and the values in each cell should be between **0 and 1** e.g.

| gene 			  | Bladder.bladder_cell  | ... | Trachea.mesenchymal_cell |
|-----------------|-----------------------|-----|--------------------------|
| ENSG00000081791 | 0.43                  | ... | 0.11                     |
| ...             | ...                   | ... | ...                      |
| ENSG00000101440 | 0.21                  | ... | 0.89                     |

During the LD score regression, the list of genes in the index column (with each gene assigned a value of 1) will be used as an indepedent variable.

In the `config.yml` pleas provide this as what you would like your output to be named and a path to the matrix.

###GWAS summary statistics

This must be one of several already munged (using the `munge-sumstats.py` script found in LD score regression) summary statistics for a given trait.

###Output directory

This is a path to a directory where you would like your output to be saved. Ideally use a path to a SSD to speed up computation. **2-3 GBs of space are usually needed** for each Specificity Input but additional GWAS summary stats will not take up much more storage.



When the config file contains the above we can run CELLECT-LDSC by navigating to the cloned CELLECT directory and entering the following:

```bash
snakemake --use-conda
```
If you wish to use multiple threads you can pass the `-j` option followed by a number to specify how many.

### Example

The config file is pre-configured to run an example which uses two expression specificity inputs created with [CELLEX](https://github.com/perslab/CELLEX) from two single cell expression atlases: [Mouse Nervous System](https://www.sciencedirect.com/science/article/pii/S009286741830789X) and Tabula Muris(https://www.nature.com/articles/s41586-018-0590-4).

The GWAS input are the UK Biobank's 1.1 million individual study on [Educational Attainment](https://www.nature.com/articles/s41588-018-0147-3) and 674 thousand individual study on Height (refrence missing).

If you have downloaded the CELLECT data then the config file is already set-up to run this example - it should take about an hour or two to complete.
