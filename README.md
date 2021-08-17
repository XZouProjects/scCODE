# scCODE
scCODE (Consensus Optimization of Differentially Expressed gene detection for single cell)
## Installation
#Package preparation
install.packages("scCODE_0.0.0.9000.tar.gz", repos = NULL, type="source")
### Cran packages
necessary1 <- c('doParallel', 'samr','doSNOW','pls','pheatmap')

installed <- necessary1 %in% installed.packages()[, 'Package']

if (length(necessary1[!installed]) >=1){
  install.packages(necessary1[!installed])
}
### Bioconductor packages
necessary2<-c('DESeq2', 'DEsingle', 
              'edgeR', 'limma', 'MAST', 'S4Vectors', 'scDD', 'scmap', 'SingleCellExperiment', 'SummarizedExperiment')
              
installed <- necessary2 %in% installed.packages()[, 'Package']

if (length(necessary2[!installed]) >=1){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install(necessary2[!installed])
}
### Github package BPSC &OGFSC，download on the site below
#### OGFSC download @ https://github.com/XZouProjects/OGFSC-R/blob/master/OGFSC_0.2.3.tar.gz
#### BPSC @ https://github.com/nghiavtr/BPSC/releases/tag/v0.99.2

install.packages("BPSC_0.99.2.tar.gz", repos = NULL, type="source")

install.packages("OGFSC_0.2.3.tar.gz", repos = NULL, type="source")

## Install scCODE

install.packages("scCODE_1.0.0.0.tar.gz", repos = NULL, type="source")

## Run scCODE

### load package

library(scCODE)

### load data
##This should be a matrix n genes by N cells, rowname should be gene names, colnames shall be cell names

datarow<-data_sample 

###The group information of cells, a factor. 1 for celltype 1, 2 for cell type 2 
group<-group_sample

### run CODE
###light = True, run CODE in a light version, which saves time.

results<-run_CODE(data_sample,group_sample,light = TRUE)

###The results is a list, which contains the AUCC and CDO matrix of methods, a conclusion of the optimal filtering method and DE method, and the DE results(DE gene name, p-adjust,fc) of the optimal method.

### Plot CODE
###Plot results, evaluation heatmap

CODE_plot(results)

## Other information
### Type ?run_CODE for more information

?run_CODE

### You can also use any component function (DE function or filter function) individually

?CODE.BPSC (t_test,MAST…) ###for help of DE methods

?CODE.filter_OGFSC  ### for help of filtering methods.
