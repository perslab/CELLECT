# CELLECT

**CELL**-type **E**xpression-specific integration for **C**omplex **T**raits (**CELLECT**) is a computational toolkit for  identifing likely etiologic cell-types underlying complex traits. CELLECT leverages existing genetic prioritization models to integrate single-cell transcriptomic and human genetic data when identifing likely etiologic cell-types. 

<!--- v1
![fig-integration](https://user-images.githubusercontent.com/5487016/62281981-0cb33d00-b44f-11e9-8c0b-24aaa2b7d286.png)
--->

<p align="center">
    <img src="https://user-images.githubusercontent.com/5487016/72677599-80e5a980-3a9e-11ea-95c4-ac645cb364a4.png" width="600"/>
</p>



## How does CELLECT work?

CELLECT quantifies the association between common polygenetic GWAS signal (heritability) and cell-type expression specificity (ES) of genes using established genetic prioritization models such as [S-LDSC](https://github.com/bulik/ldsc) ([Finucane et al., 2015](https://www.nature.com/articles/ng.3404)) and [MAGMA](http://ctglab.nl/software/magma) covariate analysis ([Skene et al., 2018](https://www.nature.com/articles/s41588-018-0129-5)). The output of CELLECT is a list of prioritized etiologic cell-types for a given human complex disease or trait.

![fig-CELLECT-regression](https://user-images.githubusercontent.com/5487016/72679543-919f1b00-3ab0-11ea-8d1d-f756fe46a4a0.png)

CELLECT takes as input GWAS data and cell-type expression specificity estimates. In order to compute robust estimates of ES, we developed the computational method called **[CELLEX](https://github.com/perslab/CELLEX)** (**CELL**-type **EX**pression-specificity). CELLEX is built on the observation that different ES metrics provide complementary cell-type expression specific profiles. Our method incorporates a ‘wisdom of the crowd’ approach by integrating multiple ES metrics to obtain improved robustness and a more expressive ES measure that captures multiple aspects of expression specificity.  

<p align="center">
    <img src="https://user-images.githubusercontent.com/5487016/72679609-299d0480-3ab1-11ea-8b05-5c1678ec270a.png" width="600"/>
</p>



*Figure legend: conceptual illustration of CELLECT and CELLEX. The bottom layer shows a disease or trait with multiple genetic components (G1-G4). CELLECT integrates disease heritability estimates with cell-type expression specificity to identify the etiologic cell-types (T1 and T4) underlying the genetic components (G1 and G4). CELLEX estimates expression specificity from single-cell transcriptomic data.*

## Update log

We have implemented CELLECT-LDSC that uses [LDSC](https://github.com/bulik/ldsc) for genetic prioritization. We expect to release CELLECT-[MAGMA](https://ctg.cncr.nl/software/magma) in the near future. We also expect to develop CELLECT-[RolyPoly](https://github.com/dcalderon/rolypoly) and CELLECT-[DEPICT](https://data.broadinstitute.org/mpg/depict/). See the [CHANGELOG](https://github.com/perslab/CELLECT/blob/master/CHANGELOG.md) for details.

## Installation

**Step 1: Install git lfs**  
We use [`git lfs`](https://git-lfs.github.com/) to store the [CELLECT data files](https://github.com/perslab/CELLECT/data) on github. To download the files you need to have `git lfs` setup before you clone the repository.
On OSX: `brew install git-lfs; git lfs install` or Ubuntu:`sudo apt-get install git-lfs; git lfs install`. For other operating systems, follow [this guide](https://github.com/git-lfs/git-lfs/wiki/Installation).

**Step 2: Clone CELLECT repository**  
```
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```
The `--recurse-submodules` is needed to clone the [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) 'ldsc' ([pascaltimshel/ldsc](https://github.com/pascaltimshel/ldsc)), which is a modfied version of the original ldsc repository.
(Cloning the repo might take few minutes as the CELLECT data files (> 1-3 GB) will be downloaded. To skip downloading the data files, use `GIT_LFS_SKIP_SMUDGE=1 git clone --recurse-submodules https://github.com/perslab/CELLECT.git` instead.)

**Step 3: Install Snakemake via conda**  

CELLECT uses the workflow management software [**Snakemake**](https://snakemake.readthedocs.io/en/stable/). To make things easier for you, CELLECT snakemake workflow utilises **conda environments** to avoid any issues with software dependencies and versioning. CELLECT snakemake workflow will automatically install all necessary dependencies. All you need to do is to [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) (if conda is not already present on your system) and then install snakemake:
   
```bash
conda install -c bioconda -c conda-forge snakemake=5.4.5
```

## Getting started with CELLECT-LDSC

An example configuration file is provided, but requires additional downloads and pre-processing. In order to run the example, please follow the [CELLECT LDSC Tutorial](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial).


1. **Modify the `config-ldsc.yml` file**: specify the input GWAS summary stats and CELLEX cell-type expression specificity. These must be in the correct format - see tutorial for example.

2. **Run CELLECT-LDSC workflow**:

```bash
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config-ldsc.yml
```
We recommend running with `-j` as it will use all available cores. Specifying `-j 4` will use up to 4 cores. 

3. **Inspect the output**:

```<BASE_OUTPUT_DIR>/results/prioritization.csv```
gives you cell-type prioritization results. You can plot the .csv file to make similar plots to this:

![ALT_TEXT](https://github.com/perslab/CELLECT/blob/master/misc/CELLECT_BMI_Tabula_Muris.gif)


### CELLECT-LDSC Tutorial: 

See our [**Github wiki**](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial) for the CELLECT-LDSC tutorial.

## Getting started with CELLECT-MAGMA

COMMING SOON!

## Documentation

Please see our [**Github wiki**](https://github.com/perslab/CELLECT/wiki) for full documentation of CELLECT. The Supplementary Notes in [Timshel (bioRxiv, 2020): _Mapping heritability of obesity by cell types_](https://www.biorxiv.org/content/10.1101/2020.01.27.920033v1) also contains relevant information on the methodology.

## Acknowledgements

We gratefully acknowledge the developers of the genetic prioritization tools used in  CELLECT: [LDSC](https://github.com/bulik/ldsc) and [MAGMA](http://ctglab.nl/software/magma). In particular, Christiaan de Leeuw and Steven Gazal for their generous support. 


## Authors

- Pascal Nordgren Timshel (University of Copenhagen) [@ptimshel](https://twitter.com/ptimshel)
- Tobi Alegbe (University of Copenhagen/University of Cambridge) [@tobialegbe](https://twitter.com/tobialegbe)
- Ben Nielsen (University of Copenhagen)

## Contact

Please create an issue on the github repo if you encounter any problems using CELLECT. 
Alternatively, you may write an email to timshel(at)sund.ku.dk

## Reference
If you find CELLECT useful for your research, please consider citing the paper:  
**[Timshel (bioRxiv, 2020): _Mapping heritability of obesity by cell types_](https://www.biorxiv.org/content/10.1101/2020.01.27.920033v1)**

