# CELLECT

**CELL**-type **E**xpression-specific integration for **C**omplex **T**raits (**CELLECT**) is a computational toolkit for identifying the probability of a cell type being involved in the genetic portion of complex traits. CELLECT uses a variety of models to assign genetic signal from a genome wide association study (GWAS) to single cell transcriptomic data. To integrate the two types of data, we must first calculate how important or specific each gene is to each cell type - the Expression Specificity (ES) of the cell type. Our group has developed a method for doing this called **CELL**-type **EX**pression specificity (**CELLEX**) which can be found [here](https://github.com/perslab/CELLEX).


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
