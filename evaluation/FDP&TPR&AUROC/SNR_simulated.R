SNR_gamma=function(rate,dep){
  # rate=rate
  qrdata=read.csv('./qrdata.csv')
  
  # wilcox_res=NULL
  # ttest_res=NULL
  # BPSC_res=NULL
  # MAST_res=NULL
  # # monocle_res=NULL
  # limma_res=NULL
  # edgeR_res=NULL
  # DEsingle_res=NULL
  # samr_res=NULL
  # DESeq2_res=NULL
  # scDD_res=NULL
  for (i in 1:length(qrdata$Data)) {
    #i=1
    snr_mean=NULL
    snr_cv=NULL
    snr_fz=NULL
    snr_var=NULL
    id_de=0
    
    celltypenames<-c(qrdata$cell1[i],qrdata$cell2[i])
    dataset<-paste0(qrdata$Data[i],'.rds')
    typ<-qrdata$typ[i]
    typeforFPR<-qrdata$tf[i]###one cell type used for null datasets building
    nbin=qrdata$nbins[i]
    
    suppressPackageStartupMessages(library(SummarizedExperiment))
    suppressPackageStartupMessages(library(MultiAssayExperiment))
    library(dplyr)
    library(ggfortify)
    dataset_name<-substr(dataset,1,nchar(dataset)-4)
    
    (dataset <- readRDS(paste0('C:/R project/DEmethods/data/',dataset)))
    experiments(dataset)
    anota<-colData(dataset)
    (dataset_gene <- experiments(dataset)[["gene"]])
    if(typ==2){
      # sam<-anota$source_name_ch1
      sam<-anota$cell_cycle_stage  #type II
    }else{sam<-anota$source_name_ch1}
    if(dataset_name=='GSE102299-smartseq2'){
      sam=anota$treatment_protocol_ch1
    }
    if(dataset_name=='GSE94383'){
      sam=anota$characteristics_ch1.2
    }
    if(dataset_name=='EMTAB3929'){
      sam=anota$Characteristics.treatment.
    }
    if(dataset_name=='SRP073808'){
      sam=anota$LibraryName
    }
    if(dataset_name=='GSE78779-GPL17021'){
      cop=anota$characteristics_ch1.4
      typeforFPR=levels(cop)[3]
    }else{
      if(dataset_name=='GSE60749-GPL13112'|dataset_name=='GSE100911'|dataset_name=='GSE94383'){
        cop<-anota$characteristics_ch1.1
      }else{
        cop<-anota$characteristics_ch1
      }
    }
    if(dataset_name=='EMTAB3929'){
      cop=anota$Characteristics.developmental.stage.
    }
    
    dataset_gene_lstm=assays(dataset_gene)[['count']]
    
    
    #dataset_gene_lstm<-assays(dataset_gene)[['TPM']]
    data_matrix<-dataset_gene_lstm[,sam==celltypenames[1]|sam==celltypenames[2]]
    cop<-cop[sam==celltypenames[1]|sam==celltypenames[2]]
    sam<-sam[sam==celltypenames[1]|sam==celltypenames[2]]
    sname<-colnames(data_matrix)
    #all zero gene kick out
    data_matrix<-subset(data_matrix,rowSums(data_matrix)>0)
    if(typ==3){
      colnames(data_matrix)<-cop
      datarow<-data_matrix[,cop==typeforFPR]
    }else{
      
      colnames(data_matrix)<-sam
      datarow<-data_matrix[,sam==typeforFPR]# data for FPR
    }
    sname<-sname[sam==typeforFPR]##for checking error
    for (l in c(1:4)) {
      # l=2
      
      filter_method=l
      
      auroc_res=NULL
      
      
      #only contains the two selected cell group data
      warnDef<-options("warn")$warn
      warnRead<-options(warn = -1)
      library(snow)
      library(pls)
      library(OGFSC)
      #instances building & analysing
      ncells=dim(datarow)[2]
      maxsz<-floor(ncells*0.5)
      
      #maxnum group 1 instance & analyse 
      
      set.seed(123)
      
      maxlis<-sample(1:ncells,2*maxsz,replace = F)
      
      instance_maxsz<-datarow[,maxlis]
      
      ####gene num in raw data
      tot0=dim(instance_maxsz)[1]
      
      dep=dep
      des=rep(0,tot0)
      degeneidx=sample(tot0,floor(tot0*dep))
      des[degeneidx]=1
      
      
      # rate=3
      set.seed(123)
      smp=sample(100,length(degeneidx),replace = T)
      
      fc=rgamma(length(smp),shape = 4,rate = rate)
      fc[which(smp<=50)]=1/fc[which(smp<=50)]
      
      
      
      group_maxsz<-factor(c(rep('1',maxsz),rep('2',maxsz)))
      idx2=which(group_maxsz==group_maxsz[1])
      
      instance_maxsz[degeneidx,idx2]=fc*instance_maxsz[degeneidx,idx2]
      
      datarow=instance_maxsz
      
      if(filter_method==1){
        logmax<-log2(datarow+1)
        OGF = OGFSC(logmax,nBins = nbin, paral_option = 1, plot_option = 1,alpha=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999))
        idx = OGF$OGFSC_idx #gene selected after filtering by OGFSC
      }else{
        if(filter_method==2){
          if(dataset_name=='GSE81076-GPL18573'|dataset_name=='GSE62270-GPL17021'){
            idx<-which(rowSums(datarow>0)/dim(instance_maxsz)[2]>0.05)}else{
              idx=which(rowSums(datarow>0)/dim(datarow)[2]>0.25)
            } 
        }else if(filter_method==3){idx=which(rowSums(datarow)>0)}
        else{
          source('C:/R project/DEmethods/scmap_filter.R')
          idx=scmap_filter(datarow = datarow)
        }}
      #DE analyse
      datarow=datarow[idx,]
      des=des[idx]
      
      
      
      idxzzz=which(colSums(datarow)==0)
      if(length(idxzzz!=0)){
        datarow=datarow[,-idxzzz]
        group_maxsz=group_maxsz[-idxzzz]
        
      }
      ncells=dim(datarow)[2]
      idx2=which(group_maxsz==unique(group_maxsz)[2])
      logdata=log2(datarow+1)
      rnames=rownames(datarow)
      colnames(datarow)=paste0('c',seq(dim(datarow)[2]))
      
      regu=function(x){
        tx<-x[idx2]
        ty<-x[-idx2]
        temp=mean(tx)-mean(ty)
        if(temp>=0){
          temp=1
        }else{
          temp=-1
        }
        return(temp)
      }
      regulate=apply(datarow, 1, regu)
      total_gene_num<-dim(datarow)[1]
      ########
      ##111111
      #ttest
      # 
      ttestfun<-function(x){
        tx<-x[idx2]
        ty<-x[-idx2]
        result<-t.test(tx,ty)
        pvalue<-result$p.value
        return(pvalue)
      }
      
      
      results <- apply(logdata, 1, ttestfun)
      ###results=p.adjust(results)   null data sets type I error control
  
        denum=min(length(which(results<=.05)),500)
        results=p.adjust(results,method = 'fdr')
        idx1=which(results<=.05)
        #SNR
        id_de=id_de+1
        m_m=apply(datarow,1,logmean)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_mean[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,logvar)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_var[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,cv)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_cv[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,fz)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_fz[id_de]=snr(m_s,m_ns)
        #write.csv(DE_results1,paste0(filename,rate,'FC-ttest-adj_results.csv'))
   
      # # #######
      # ##222222
      #   # ####BPSC
      library(BPSC)
      cl=makeCluster(18)
      registerDoParallel(cl)
      
      
      controlIds=which(group_maxsz==1)
      design=model.matrix(~group_maxsz)
      coef=2
      res=BPglm(data=datarow, controlIds=controlIds, design=design, coef=coef, estIntPar=F, useParallel=T)
      ###results=p.adjust(results)   null data sets type I error control
      total_gene_num<-dim(datarow)[1]
        results=res$PVAL
        denum=min(length(which(results<=.05)),500)
        results=p.adjust(results,method = 'fdr')
        idx1=which(results<=.05)
        #DE_p<-length(idx)/total_gene_num
        #SNR
        id_de=id_de+1
        m_m=apply(datarow,1,logmean)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_mean[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,logvar)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_var[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,cv)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_cv[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,fz)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_fz[id_de]=snr(m_s,m_ns)
        #write.csv(DE_results1,paste0(filename,rate,'FC-BPSC-adj_results.csv'))
        stopCluster(cl)
  
      
      # # ###
      # # ##33333
      # ##DESeq2
      
      library(DESeq2)
      datarow_in<-apply(datarow, 2, as.integer)
      
      dds_ins_maxsz<-DESeqDataSetFromMatrix(datarow_in,DataFrame(group_maxsz),~group_maxsz)
      des_ins_maxsz<-try(DESeq(dds_ins_maxsz,parallel = F),silent = T)
      if('try-error' %in% class(des_ins_maxsz)){
        fakegene=round(colSums(datarow_in)/dim(datarow_in)[1])+1
        datarow_in=rbind(datarow_in,fakegene)
        dds_ins_maxsz<-DESeqDataSetFromMatrix(datarow_in,DataFrame(group_maxsz),~group_maxsz)
        des_ins_maxsz<-try(DESeq(dds_ins_maxsz,parallel = F),silent = T)
      }
      res_ins_maxsz<-DESeq2::results(des_ins_maxsz)
      
      
      ###results=p.adjust(results)   null data sets type I error control
      
        results=res_ins_maxsz$pvalue
        denum=min(length(which(results<=.05)),500)
        results=p.adjust(results,method = 'fdr')
        if(length(results)>length(des)){
          results=results[-length(results)]
        }
        idx1=which(results<=.05)
        #DE_p<-length(idx)/total_gene_num
        #SNR
        id_de=id_de+1
        m_m=apply(datarow,1,logmean)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_mean[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,logvar)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_var[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,cv)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_cv[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,fz)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_fz[id_de]=snr(m_s,m_ns)
        #write.csv(DE_results1,paste0(filename,rate,'FC-DESeq2-adj_results.csv'))
      
      
      # ###MAST
      # ##44444
      # #
      library(MAST)
      cData<-data.frame(wellKey=paste0('C',1:dim(logdata)[2]))
      fData<-data.frame(primerid=rownames(logdata))
      colnames(logdata)=paste0('C',1:dim(logdata)[2])
      sca<-FromMatrix(logdata,cData,fData)
      colData(sca)$cond<-group_maxsz
      zlmcond<-try(zlm(~cond,sca,parallel = T),silent = T)
      summarycond<-summary(zlmcond,doLRT='cond2')
      summaryDt <- summarycond$datatable
      fcHurdle <- merge(summaryDt[contrast=='cond2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                        summaryDt[contrast=='cond2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
      
      ###results=p.adjust(results)   null data sets type I error control
      
        results=fcHurdle$`Pr(>Chisq)`
        denum=min(length(which(results<=.05)),500)
        results=p.adjust(results,method = 'fdr')
        idx1=which(results<=.05)
        #DE_p<-length(idx)/total_gene_num
        #SNR
        id_de=id_de+1
        m_m=apply(datarow,1,logmean)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_mean[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,logvar)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_var[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,cv)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_cv[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,fz)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_fz[id_de]=snr(m_s,m_ns)
        #write.csv(DE_results1,paste0(filename,rate,'FC-MAST-adj_results.csv'))
      
      
      
      # ###wilcox test
      # ###555555
      # # ##
      wilcoxtestfun<-function(x){
        tx<-x[idx2]
        ty<-x[-idx2]
        result<-wilcox.test(tx,ty)
        pvalue<-result$p.value
        return(pvalue)
      }
      
      
      results <- apply(logdata, 1, wilcoxtestfun)
      ###results=p.adjust(results)   null data sets type I error control
      total_gene_num<-dim(datarow)[1]
    
        
        denum=min(length(which(results<=.05)),500)
        results=p.adjust(results,method = 'fdr')
        idx1=which(results<=.05)
        #DE_p<-length(idx)/total_gene_num
        #SNR
        id_de=id_de+1
        m_m=apply(datarow,1,logmean)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_mean[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,logvar)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_var[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,cv)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_cv[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,fz)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_fz[id_de]=snr(m_s,m_ns)
        #write.csv(DE_results1,paste0(filename,rate,'FC-wilcoxtest-adj_results.csv'))
      
      
      # ###DEsingle
      # ##66666
      # ##
      library(DEsingle)
      results <- try(DEsingle(counts = datarow,group = group_maxsz,parallel = T),silent = T)
      if('try-error'%in% class(results)){
        id_de=id_de+1
        snr_mean[id_de]=NA
        snr_var[id_de]=NA
        snr_cv[id_de]=NA
        snr_fz[id_de]=NA
        #write.csv(DE_results1,paste0(filename,rate,'FC-DEsingle-adj_results.csv'))
      }else{
        results=results$pvalue
        
        denum=min(length(which(results<=.05)),500)
        results=p.adjust(results,method = 'fdr')
        idx1=which(results<=.05)
        #DE_p<-length(idx)/total_gene_num
        #SNR
        id_de=id_de+1
        m_m=apply(datarow,1,logmean)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_mean[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,logvar)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_var[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,cv)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_cv[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,fz)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_fz[id_de]=snr(m_s,m_ns)
        #write.csv(DE_results1,paste0(filename,rate,'FC-DEsingle-adj_results.csv'))
        
      }
      #   #
      #   # ###edgeR
      #   # ##7777777
      #   # #
      cl=makeCluster(16)
      registerDoParallel(cl)
      library(edgeR)
      
      dgelist =DGEList(counts = datarow, group = group_maxsz)
      dgelist_norm <- try(calcNormFactors(dgelist, method = 'TMM'),silent=T)
      if('try-error'%in%class(dgelist_norm)){
        id_de=id_de+1
        snr_mean[id_de]=NA
        snr_var[id_de]=NA
        snr_cv[id_de]=NA
        snr_fz[id_de]=NA
        #write.csv(DE_results1,paste0(filename,rate,'FC-edgeR-adj_results.csv'))
      }else{
        design=model.matrix(~group_maxsz)
        dge <- estimateDisp(dgelist_norm, design, robust = F) #dispersion not robust
        fit <- glmFit(dge, design, robust = F)     #fit
        lrt <- glmLRT(fit)   #test
        results=lrt$table$PValue
        
        denum=min(length(which(results<=.05)),500)
        results=p.adjust(results,method = 'fdr')
        idx1=which(results<=.05)
        #DE_p<-length(idx)/total_gene_num
        #SNR
        id_de=id_de+1
        m_m=apply(datarow,1,logmean)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_mean[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,logvar)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_var[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,cv)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_cv[id_de]=snr(m_s,m_ns)
        m_m=apply(datarow,1,fz)
        m_s=m_m[idx1]
        m_ns=m_m[-idx1]
        snr_fz[id_de]=snr(m_s,m_ns)
        #write.csv(DE_results1,paste0(filename,rate,'FC-edgeR-adj_results.csv'))
      }
      stopCluster(cl)
      #
      
      # #samr
      library(samr)
      data=list(x=datarow,y=group_maxsz,geneid=as.character(1:nrow(datarow)),
                genenames=rownames(datarow),logged2=TRUE)
      res=samr(data,resp.type = 'Two class unpaired',nperms = 1000)
      results=samr.pvalues.from.perms(res$tt, res$ttstar)
      denum=min(length(which(results<=.05)),500)
      results=p.adjust(results,method = 'fdr')
      idx1=which(results<=.05)
      #DE_p<-length(idx)/total_gene_num
      #SNR
      id_de=id_de+1
      m_m=apply(datarow,1,logmean)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_mean[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,logvar)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_var[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,cv)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_cv[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,fz)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_fz[id_de]=snr(m_s,m_ns)
      #write.csv(DE_results1,paste0(filename,rate,'FC-samr-adj_results.csv'))
      
      # ###
      # ##10101010
      # #limma
      library(limma)
      limmagroup<-c(rep('control',dim(datarow)[2]))
      limmagroup[idx2]='fake'
      limmagroup=factor(limmagroup)
      names(limmagroup)=colnames(datarow)
      design=model.matrix(~0+limmagroup)
      colnames(design)=levels(factor(limmagroup))
      rownames(design)=colnames(datarow)
      norm=voom(datarow,design =design,plot = FALSE)
      fit=lmFit(norm,design,method = 'ls')
      contrast=makeContrasts('control-fake',levels = design)
      fit2=contrasts.fit(fit,contrast)
      fit2=eBayes(fit2)
      results=fit2$p.value
      
      denum=min(length(which(results<=.05)),500)
      results=p.adjust(results,method = 'fdr')
      idx1=which(results<=.05)
      #DE_p<-length(idx)/total_gene_num
      #SNR
      id_de=id_de+1
      m_m=apply(datarow,1,logmean)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_mean[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,logvar)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_var[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,cv)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_cv[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,fz)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_fz[id_de]=snr(m_s,m_ns)
      #write.csv(DE_results1,paste0(filename,rate,'FC-limma-adj_results.csv'))
      
      #   ####scDD
      library(scDD)
      library(SingleCellExperiment)
      names(group_maxsz)=colnames(datarow)
      condition=group_maxsz
      idx=which(colSums(datarow)==0)
      if(length(idx)!=0){
        datarow=datarow[,-idx]
        condition=condition[-idx]
      }
      
      sce=SingleCellExperiment(assays=list(counts=datarow),colData=data.frame(condition))
      
      scDatEx.scran <- preprocess(sce, zero.thresh=1.0, scran_norm=TRUE)
      
      prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
      scDatExSim <- scDD(scDatEx.scran, prior_param=prior_param, categorize = F,testZeroes=FALSE)
      
      res_scDD=scDD::results(scDatExSim)
      results=res_scDD$nonzero.pvalue
      
      results=p.adjust(results,method = 'fdr')
      idx1=which(results<=.05)
      #DE_p<-length(idx)/total_gene_num
      #SNR
      id_de=id_de+1
      m_m=apply(datarow,1,logmean)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_mean[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,logvar)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_var[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,cv)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_cv[id_de]=snr(m_s,m_ns)
      m_m=apply(datarow,1,fz)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_fz[id_de]=snr(m_s,m_ns)
    }
    #
    #
    SNR_res=rbind(snr_mean,snr_var,snr_cv,snr_fz)
    rownames(SNR_res)=c('mean','var','cv','fz')
    colnames(SNR_res)=rep(c('ttest','BPSC','DESeq2','MAST','wilcoxtest','DEsingle','edgeR','samr','limma','scDD'),4)
    write.csv(SNR_res,paste0('C:/R project/DEmethods/CODE/simulation/',qrdata$Data[i],'/SNR_res.csv'))
    
    # colnames(monocle_res)=c('DE_p','fc','totalgene','total_gene_afterfilter','trueDEnum_after_filter','DEgenenum','fdr','tpr','auroc')
    # write.csv(monocle_res,paste0(filename,rate,'FC-monocle-results.csv'))
    cat('complete analysing dataset x!\n')
    cat('\n ncells',ncells)
  }
}
