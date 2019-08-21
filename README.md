# CELLECT

**CELL**-type **E**xpression-specific integration for **C**omplex **T**raits (**CELLECT**) is a computational toolkit for  identifing likely etiologic cell-types underlying complex traits. CELLECT leverages existing genetic prioritization models to integrate single-cell transcriptomic and human genetic data when identifing likely etiologic cell-types. Currently we have only implemented CELLECT-LDSC that uses [LDSC](https://github.com/bulik/ldsc) for genetic prioritization. We expect to release CELLECT-MAGMA, CELLECT-RolyPoly and CELLECT-DEPICT in the near future.

<!---
![fig-integration](https://user-images.githubusercontent.com/5487016/62281981-0cb33d00-b44f-11e9-8c0b-24aaa2b7d286.png)
--->

<p align="center">
    <img src="https://user-images.githubusercontent.com/5487016/62281981-0cb33d00-b44f-11e9-8c0b-24aaa2b7d286.png" width="600"/>
</p>



## How does CELLECT work?

CELLECT quantifies the association between common polygenetic GWAS signal (heritability) and cell-type expression specificity (ES) of genes using established genetic prioritization models such as LDSC (Hilary Kiyo Finucane et al., 2015), RolyPoly (Calderon et al., 2017), DEPICT (Pers et al., 2015) or MAGMA covariate analysis (Skene et al., 2018).

CELLECT takes as input GWAS data and cell-type expression specificity estimates. In order to compute robust estimates of ES, we developed the computational method called **CELLEX** (**CELL**-type **EX**pression-specificity). CELLEX is built on the observation that different ES metrics provide complementary cell-type expression specific profiles. Our method incorporates a ‘wisdom of the crowd’ approach by integrating multiple ES metrics to obtain improved robustness and a more expressive ES measure that captures multiple aspects of expression specificity.  CELLEX can be found [here](https://github.com/perslab/CELLEX).

![fig-CELLECT-conceptual-h2](https://user-images.githubusercontent.com/5487016/62367093-e3ff7600-b528-11e9-8879-8f69005fbea5.png)

Schematic illustration of CELLECT and CELLEX. The bottom layer shows a disease or trait with multiple genetic components (G1-G4). CELLECT integrates disease heritability estimates with cell-type expression specificity to identify the etiologic cell-types (T1 and T4) underlying the genetic components (G1 and G4). CELLEX estimates expression specificity from single-cell transcriptomic data.


## Installation

1. **Install git lsf**
We use [`git lfs`](https://git-lfs.github.com/) to store the [CELLECT data files](https://github.com/perslab/CELLECT/data) on github. To download the files you need to have `git lfs` setup before you clone the repository.
Mac:
```
brew update
brew install git-lfs
git lfs install
```
Ubuntu:
```
sudo apt-get install git-lfs
git lfs install
```
For other operating systems, follow [this guide](https://github.com/git-lfs/git-lfs/wiki/Installation).


2. **Clone CELLECT repository**
```
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```
The `--recurse-submodules` is needed to clone the [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) 'ldsc' ([pascaltimshel/ldsc](https://github.com/pascaltimshel/ldsc)), which is a modfied version of the original ldsc repository.
(This might take few minutes as the CELLECT data files (> 1 GB) will be downloaded using git lfs. To skip downloading the data files, use `GIT_LFS_SKIP_SMUDGE=1 git clone --recurse-submodules https://github.com/perslab/CELLECT.git` instead.)

3. **Install Snakemake via conda**

CELLECT uses the workflow management software [**Snakemake**](https://snakemake.readthedocs.io/en/stable/). To make things easier for you, CELLECT snakemake workflow utilises **conda environments** to avoid any issues with software dependencies and versioning. CELLECT snakemake workflow will automatically install all necessary dependencies. All you need to do is to [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) (if conda is not already present on your system) and then install snakemake:
   
```bash
conda install -c bioconda snakemake
```

## Getting started with CELLECT-LDSC


1. **Modify the `config-ldsc.yml` file**. The config file is divided up into two categories:

* **RUN-SPECIFIC**: These variables can change from run-to-run and affect how CELLECT processes its input data.
* **CONSTANTS**: These variables include things like paths to CELLECT data and libraries that generally do not need to be changed. 


2. When the config file contains the above we can run CELLECT-LDSC by navigating to the cloned CELLECT directory and entering the following:

```bash
snakemake --use-conda -j 4 -s cellect-ldsc.snakefile
```
This will run the workflow using 4 cores (`-j 4`). If you wish to use to use all available cores pass just the `-j` flag.


### CELLECT-LDSC Example: 

This example with run CELLECT on two 
The example should take 5-15 minutes to run depending on the available number of cores on your system.
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

Please see our [**github wiki**](https://github.com/perslab/CELLECT/wiki) for additional documentation.

## Authors

Pascal Nordgren Timshel (University of Copenhagen)
Tobi Alegby (University of Copenhagen)
Ben Nielsen (University of Copenhagen)


