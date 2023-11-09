![图片](https://github.com/XZouProjects/scCODE/assets/17633478/fc6837b4-5e64-48eb-9ff3-2811711e609b)![图片](https://user-images.githubusercontent.com/17633478/137343572-3b77beaf-d70e-4001-bd6a-fe27fd3f2628.png)
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

    install.packages("scCODE_1.2.0.0.tar.gz", repos = NULL, type="source")

## Run scCODE

We can start using scCODE by running the sample data.

    library(scCODE)

The input requires a count matrix (genes by cells), data_sccode contains the CD4+ T cells (naive and activated).

    data1<-data1_sccode
    
    data2<-data2_sccode

There are two parameters selectable. "light", True or False, default as True, run scCODE in a light version which saves time. "top_ranked", the number of top-ranked strategies selected (5-10), default as 5.

Now, we can run the sample data like below:

    results<-scCODE(data1,data2,light = TRUE,top_ranked=5)

### Output results

The results is a list containing the DE gene results and the metrics evaluation results.

The DE results is a dataframe which contains the ranked DE gene name along with its times-detected, logFC (group1/group2, and group 1 refers to the unique(group)[1]). it can be access by :

    results$DE_results
    
The metrics results can be accessed by the code below, which is the Z-score of the metrics of each analysis strategies. 

    results$`Z-score-metrics`

### Lotus plot, divided by detected_times
    
    DE_results<-results$DE_results

###times by detected times
    times<- 1

    DE_results$logFC_by_D<-ifelse(DE_results$logFC>=0,DE_results$logFC+(DE_results$Detected_times-1)*times,DE_results$logFC-(DE_results$Detected_times-  1)*times)
    DE_results$Detected_times<-as.factor(DE_results$Detected_times)
    library(ggplot2)
    library(RColorBrewer)
    coul <- rev(brewer.pal(5, 'Set1'))

    P<-ggplot(
      DE_results, 
      aes(x = logFC_by_D, 
          y = -log10(P_adjust), 
          colour=Detected_times)) +
      geom_point(alpha=0.4, size=1.4) +
      scale_color_manual(values=coul)+
      labs(x="log2(fold change_by_detimes)",
           y="-log10 (p-value)")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position="right", 
            legend.title = element_blank()
      )
    
    tiff('Fig lotus_sccode.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
    P
    dev.off()


