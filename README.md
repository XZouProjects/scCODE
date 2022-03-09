![图片](https://user-images.githubusercontent.com/17633478/137343572-3b77beaf-d70e-4001-bd6a-fe27fd3f2628.png)
# scCODE

scCODE (Consensus Optimization of Differentially Expressed gene detection for single cell)

## Installation

Most dependent packages can be installed by the code below:

    necessary1 <- c('doParallel', 'samr','doSNOW','pls')

    installed <- necessary1 %in% installed.packages()[, 'Package']


    if (length(necessary1[!installed]) >=1){

      install.packages(necessary1[!installed])
  
     }

    necessary2<-c('DESeq2', 'DEsingle', 
              'edgeR', 'limma', 'MAST', 'S4Vectors', 'scDD', 'scmap', 'SingleCellExperiment', 'SummarizedExperiment')
              
              
    installed <- necessary2 %in% installed.packages()[, 'Package']


    if (length(necessary2[!installed]) >=1){

      if (!requireNamespace("BiocManager", quietly = TRUE))
  
        install.packages("BiocManager")
    
        library(BiocManager)
    
        BiocManager::install(necessary2[!installed])
    
      }
   
Two packages OGFSC and BPSC are required and need to install from Github.

First, check if we have already installed the two packages by :

    'OGFSC'%in% installed.packages()
    
    'BPSC'%in% installed.packages()
    
If not, we can download the installation file and install them through [OGFSC](https://github.com/XZouProjects/OGFSC-R/blob/master/OGFSC_0.2.3.tar.gz) and [BPSC](https://github.com/nghiavtr/BPSC/releases/tag/v0.99.2)

## Install scCODE

Now, we can install the scCODE by downloading the installation file in this page (scCODE_1.0.1.0.tar.gz), and install it:

    install.packages("scCODE_1.0.1.0.tar.gz", repos = NULL, type="source")

## Run scCODE

We can start using scCODE by running the sample data.

    library(scCODE)

The input requires a count matrix (genes by cells), data_sccode contains the CD4+ T cells (naive and activated).

    data<-data_sccode

The input also requires a vector, group, which is the cell group information of the cells. The group_sccode is the naive or activated (1 or 2) information of the data.

    group<-group_sccode


There are two parameters selectable. "light", True or False, default as True, run scCODE in a light version which saves time. "top_ranked", the number of top-ranked strategies selected (5-10), default as 5.

Now, we can run the sample data like below:

    results<-scCODE(data,group,light = TRUE,top_ranked=5)

### Output results

The results is a list containing the DE gene results and the metrics evaluation results.

The DE results is a dataframe which contains the ranked DE gene name along with its times-detected, logFC (group1/group2, and group 1 refers to the unique(group)[1]). it can be access by :

    results$DE_results
    
The metrics results can be accessed by the code below, which is the Z-score of the metrics of each analysis strategies. 

    results$Z-score-metrics
