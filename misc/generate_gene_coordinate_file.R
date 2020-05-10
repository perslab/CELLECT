############### SYNOPSIS ###################
# Create gene coordinate file compatible with CELLECT-LDSC and CELLECT-MAGMA

### DESCRIPTION
# ...

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:
# GRCh37 assembly in Ensembl: https://www.ensembl.org/info/website/tutorials/grch37.html
# https://grch37.ensembl.org/biomart 
# Ensembl GRCh37: Full Feb 2014 archive with BLAST, VEP and BioMart
# GRCh37.p13

# Ensembl biomaRt guide: http://www.ensembl.org/info/data/biomart/biomart_r_package.html


############################################

### Installation
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)
library(tidyverse)

devtools::session_info() 
# R version 3.5.3 (2019-03-11)
# biomaRt       * 2.38.0    2018-10-30 [1] Bioconductor  


# ======================================================================= #
# ==============================  biomaRt  ============================== #
# ======================================================================= #
### Set version variable
#ENSEMBL_VERSION = NULL # for newest version
ENSEMBL_VERSION = 91 # Ensembl release 91 = December 2017 | for specific version
GENOME_VERSION = 37 # GRCh37=hg19; GRCh38=hg38

### Connect to the BioMart database and dataset hosted by Ensembl
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ENSEMBL_VERSION, GRCh=GENOME_VERSION, verbose=T)

df.ensembl.attr <- listAttributes(ensembl) # ~1280 attributes

# ======================================================================= #
# ========================= PAGE = GENES FEATURES  ====================== #
# ======================================================================= #

### Get all Ensembl Gene IDs - JUST TO SEE HOW MANY GENES THERE ARE IN TOTAL
df.BM.all.ids <- getBM(attributes = c("ensembl_gene_id"), mart=ensembl)
nrow(df.BM.all.ids) # 63677 for ensembl_v91
str(df.BM.all.ids)

### Features page
df.BM.feature <- getBM(attributes = c("ensembl_gene_id", 
                                      "chromosome_name",
                                      "start_position", # Gene start (bp) [transcription start site]
                                      "end_position", # Gene end (bp) [transcription end site]
                                      "gene_biotype",
                                      "strand",
                                      "external_gene_name"), mart=ensembl)
str(df.BM.feature)
nrow(df.BM.feature) 

stopifnot(length(unique(df.BM.feature$ensembl_gene_id)) == nrow(df.BM.feature)) # --> TRUE, all ensembl_gene_id are unique


# ======================================================================= #
# ===============================  INSPECT  ============================== #
# ======================================================================= #

df.BM.feature %>% count(chromosome_name) %>% arrange(desc(n)) %>% print(n=50)
# ---> Output *does* contain 'X' and 'Y' chromosomes
df.BM.feature %>% distinct(chromosome_name) %>% nrow() # ---> 265 different chromosomes

df.BM.feature %>% distinct(strand) # --> ONLY 1 / -1

# ======================================================================= #
# ===============================  MODIFY/FILTER  ============================== #
# ======================================================================= #


### Filter
df.gene_coord <- as_tibble(df.BM.feature) %>% 
  filter(gene_biotype=='protein_coding') %>% # protein-coding only
  filter(chromosome_name %in% seq(1,22)) # keep autosomal genes only

### Modify
df.gene_coord <- df.gene_coord %>% 
  select(-gene_biotype) %>%
  mutate(strand=if_else(strand==1, "+", "-")) %>%
  arrange(chromosome_name, start_position)

df.gene_coord  

# ======================================================================= #
# ===============================  EXPORT  ============================== #
# ======================================================================= #

### Columns
# ensembl_gene_id
# chromosome
# start_position
# stop_position
# strand
# gene_name

### Write table
file.ensmbl_annotation <- "gene_coordinates.GRCh37.ensembl_v91.txt"
df.gene_coord %>% write_tsv(path=file.ensmbl_annotation, col_names=F) # no header


# ======================================================================= #
# ===============================  WIKI  ============================== #
# ======================================================================= #

# GENELOC_FILE
# Requires a formatted data file, with rows corresponding to genes; 
# a header is not allowed. 
# The file must have the four columns [in that order]: 
# gene ID
# chromosome
# start position
# stop position. 
# It may have a fifth column containing the strand, coded as + for the positive/sense strand and - for the negative/antisense strand
# Example:
# GENE_ID1   1       69091   70008   +       OR4F5
# GENE_ID2       1       142447  174392  -       LOC100996442
# GENE_ID3  1       367659  368597  +       OR4F29


# Autosomal chromosomes are coded numerically from 1 through 22, the sex chromosomes can be coded either
# as X and Y or as 23 and 24. If included, gene strands must be coded as + for the sense/positive strand
# and - for the antisense/negative strand.
