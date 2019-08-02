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


One of the models used by CELLECT is LD score regression - it is vital that our forked version of this repository [(found here)](https://github.com/pascaltimshel/ldsc) is also cloned.

The models in CELLECT are built in Snakemake and the pipelines utilise conda environments. Therefore the easiest way to get started would be to [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) (if conda is not already present on your system) and then either within base conda or a conda environment:
```bash
conda install -c bioconda -c conda-forge snakemake
conda install pandas
```

There are several data files (roughly 6GB after unpacking) which are required to use CELLECT they can be downloaded as follows:
```bash
wget www.linktocellect.com/cellect-data.gz
tar -xvf cellect-data.gz
```

Finally, within your cloned version of this repository, modify the `config.yml` file so that it is specific to the system you are working on. The config file is divided up into two categories (plus model-specific sub categories):

* **RUN-SPECIFIC**: These variables can change from run-to-run and affect how CELLECT processes its input data
* **CONSTANTS**: These variables will generally remain constant and will includes things like paths to data and scripts

For now, only modify the CONSTANTS part to point to the cloned LDSC script and the data directory we just downloaded.

## Usage of CELLECT-LDSC

To run CELLECT-LDSC we must specify in the `config.yml` either one or multiple multigeneset files as well as one or multiple GWAS summary statistics.

The multigeneset file contains Expression Specificity values for every specific gene in each cell type and can be created by providing CELLEX with single cell count matrix and metadata.

The GWAS summary statistics are pre-processed by the LDSC script `munge_sumstats.py`.

When multigeneset and summary statistic names have been saved in the config file we can run CELLECT-LDSC by navigating to the cloned CELLECT directory and entering the following:

```bash
snakemake --use-conda
```
If you wish to use multiple threads you can pass the `-j` option followed by a number to specify how many.

### Example

To demonstrate CELLECT-LDSC we have a short example that takes roughly 30 minutes to run on ... 

Our single cell input for this example is the [Mouse Nervous System](https://www.sciencedirect.com/science/article/pii/S009286741830789X) processed with CELLEX to get specificity scores then a subset of cell types taken to speed up the running time. The GWAS input is the UK Biobank's 1.1 million individual study on [Educational Attainment](https://www.nature.com/articles/s41588-018-0147-3)

If you have downloaded the CELLECT data then the config file is already set-up to run this example.
