#####splatter simulate
require(splatter)
require(scater)

setwd('/data')

###SNR function
cv=function(a){
  if(mean(a)<=0){
    res=0
  }else{
    res=var(a)/mean(a)
  }
  return(res)
}
snr=function(s,ns){
  snr_value=(mean(s)-mean(ns))/(var(s)+var(ns))
  return(snr_value)
}
fz=function(a){
  fz_value=length(which(a==0))/length(a)
}
logmean=function(a){
  value=log2(mean(a)+1)
  return(value)
}
logvar=function(a){
  value=log2(var(a)+1)
  return(value)
}

source('./New_Metrics.R')
source('./DE_functions.R')
source('./Filter_functions.R')
source('./auroc.R')

# nts<-5

  dataset_name<-'EMTAB_h_k'
  
  data_new<-readRDS(paste0(dataset_name,'.rds'))
  datarow<-data_new
  #instances building & analysing
  ncells=dim(datarow)[2]
  maxsz<-floor(ncells*0.5)
  datarow<-datarow[which(rowSums(datarow)>0),]
  # dir.create(paste0('C:/R project/DEmethods/simulate buquan/',dataset_name,'/'))
  
  group<-rep(2,ncol(datarow))
  group[which(colnames(datarow)==unique(colnames(datarow))[1])]=1
    
  group=factor(group)
    
    ### appeared later in group as cell type 2
    idx2<-which(group==unique(group)[2])
    ##cal_fc for each DE genes detected
    cal_fc<-function(data){
      ex_1<-data[-idx2]
      ex_2<-data[idx2]
      fc<-mean(ex_1)/mean(ex_2)
      return(fc)
    }
    
    
    ####gene num in raw data
    tot0=dim(datarow)[1]
    
    
    colnames(datarow)<-paste0('cell',seq(ncol(datarow)))
    
    data_ori<-datarow
    data_ori<-as.matrix(data_ori)
    # data_ori<-data_ori[which(rowSums(data_ori)>0),]
    group_ori<-group
    # nts<-stats::median(c(nts,5,10))### top 5 - 10 consensus
    
    genename_all<-rownames(data_ori)
    
    ###FC was recorded as the value of origin data input
    fc_all<-apply(data_ori, 1, cal_fc)
    
    ###log2fc
    log2fc_all<-abs(log2(fc_all))
    
    ###SNR and bias control
    logmean_all=apply(data_ori,1,logmean)
    
    logvar_all=apply(data_ori,1,logvar)
    
    cv_all=apply(data_ori,1,cv)
    
    fz_all=apply(data_ori,1,fz)
    
    
    
    
    
    DEmethods<-rev(c('scDD','MAST','t_test','BPSC','wilcox_test','samr','DESeq2','limma','DEsingle','edgeR'))
    
    filtermethods<-c(1:4)

    ######SNR record
    snr_mean_record<-NULL
    
    snr_fz_record<-NULL
    
    snr_var_record<-NULL
    
    snr_cv_record<-NULL
    
    ######mean fz var cv fc
    logmean_record<-NULL
    
    fz_record<-NULL
    
    logvar_record<-NULL
    
    cv_record<-NULL
    
    logfc_record<-NULL
    
    for(filt in filtermethods){
      if(filt==1){
        ###filtering by OGFSC
        idx<-scCODE.filter_OGFSC(data_ori)
      }else{
        if(filt==2){
          ###filtering by conquer
          idx<-scCODE.filter_conquer(data_ori)
        }else{
          if(filt==3){
            ###filtering by scmap
            idx<-scCODE.filter_scmap(data_ori)
          }else{
            ###no filter
            idx=which(rowSums(data_ori)>0)
          }
        }
      }
      
      datatest<-data_ori[idx,]
      idx_c=which(colSums(datatest)==0)
      if(length(idx_c)!=0){
        datatest=datatest[,-idx_c]
        group=group_ori[-idx_c]
      }else{
        group=group_ori
      }
      datatest<-as.matrix(datatest)
      ###gene name
      genename<-rownames(datatest)
      
      for (de in DEmethods) {
        
        
        run_func<-paste0('scCODE.',de,'(datatest,group)')
        res_temp<-try(eval(parse(text = run_func)),silent = T)
        if('try-error'%in%class(res_temp)){
          res_temp=NA
        }
        ### adjust p value<=0.05
        res_temp<-stats::p.adjust(res_temp,method = 'fdr')
        names(res_temp)<-genename[1:length(res_temp)]
        ###record all gene order
        
        ###record all DE gene results
        idx_sig=which(res_temp<=.05)
        ###
        n_gene=length(idx_sig)

        idx_sig=order(res_temp)[1:n_gene]
        ###Padj record
        padj_temp<-res_temp[idx_sig]
        ### All DE gene
        gene_all<-genename[idx_sig]
        
        ### FC record for each result, idx_FC, match back to the origin data
        idx_fc<-match(gene_all,genename_all)
        
        fc_temp<-fc_all[idx_fc]
        
        logmean_de<-logmean_all[idx_fc]
        
        logmean_record<-c(logmean_record,mean(logmean_de))
        
        logmean_non<-logmean_all[-idx_fc]
        
        snr_mean_temp=snr(logmean_de,logmean_non)
        
        snr_mean_record<-c(snr_mean_record,snr_mean_temp)
        
        logvar_de<-logvar_all[idx_fc]
        
        logvar_record<-c(logvar_record,mean(logvar_de))
        
        logvar_non<-logvar_all[-idx_fc]
        
        snr_var_temp=snr(logvar_de,logvar_non)
        
        snr_var_record<-c(snr_var_record,snr_var_temp)
        
        cv_de<-cv_all[idx_fc]
        
        cv_record<-c(cv_record,mean(cv_de))
        
        cv_non<-cv_all[-idx_fc]
        
        snr_cv_temp=snr(cv_de,cv_non)
        
        snr_cv_record<-c(snr_cv_record,snr_cv_temp)
        
        fz_de<-fz_all[idx_fc]
        
        fz_record<-c(fz_record,mean(fz_de))
        
        fz_non<-fz_all[-idx_fc]
        
        snr_fz_temp=snr(fz_de,fz_non)
        
        snr_fz_record<-c(snr_fz_record,snr_fz_temp)
        
        logfc_temp<-log2fc_all[idx_fc]
        
        logfc_record<-c(logfc_record,mean(logfc_temp))
   
        
      }
      
    }
    #####SNR_real
   
    snr_experimental=cbind(logmean_record,logvar_record,cv_record,fz_record,logfc_record,snr_mean_record,snr_var_record,snr_cv_record,snr_fz_record)
    rownames(snr_experimental)=c(paste0('OGFSC-',DEmethods),paste0('conquer-',DEmethods),paste0('scmap-',DEmethods),paste0('without_filter-',DEmethods))
    
    dir.create(paste0('./',dataset_name))
    
    write.csv(snr_experimental,paste0('./',dataset_name,'/SNR_expdata','.csv'))
    
    cat(paste0('/nAnalysis Finished'))

  
