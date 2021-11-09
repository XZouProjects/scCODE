![图片](https://user-images.githubusercontent.com/17633478/137343572-3b77beaf-d70e-4001-bd6a-fe27fd3f2628.png)
# scCODE
scCODE (Consensus Optimization of Differentially Expressed gene detection for single cell)
## Installation

R version >=3.6.3

### Cran packages
necessary1 <- c('doParallel', 'samr','doSNOW','pls','pheatmap', 'writexl')


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
### Github package BPSC & OGFSC，cheeck if these packages were installed, and if not, they could be downloaded at the site below.
#### OGFSC download @ https://github.com/XZouProjects/OGFSC-R/blob/master/OGFSC_0.2.3.tar.gz
#### BPSC @ https://github.com/nghiavtr/BPSC/releases/tag/v0.99.2

install.packages("BPSC_0.99.2.tar.gz", repos = NULL, type="source")


install.packages("OGFSC_0.2.3.tar.gz", repos = NULL, type="source")


### Install scCODE

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

###outdir, path to save results.

###nts, number of the top suitable strategies for consensus optimization, range 5-10.


results<-scCODE(datarow,group,light = TRUE,outdir="./",nts=5)


###The results will be saved in the path file, where contains an excel of consensus DE information, 3 excels of all the DE results by conquer, OGFSC and scmap filtering.

### plot CODE
###Plot results, evaluation heatmap

scCODE_plot(results)

## Other information

### More information about scCODE

?scCODE

### Information about the optional DE and gene filtering methods

?scCODE.BPSC (t_test,MAST…) ###for help of DE methods

?scCODE.filter_OGFSC  ### for help of filtering methods.
