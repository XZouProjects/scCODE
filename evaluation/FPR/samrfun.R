samr_anlys<-function(typ,dataset,celltypenames,typeforFPR,filter_method,nbin){
  
  # celltypenames<-c('16-cell stage blastomere','Mid blastocyst cell (92-94h post-fertilization)')
  # dataset<-'GSE45719.rds'
  # 
  # typeforFPR<-celltypenames[1]
  # groupsize<-c(24,12,6)
  # sOgfsc=T
  # Ogfsc=F
  if(missing(nbin)){nbin<-20}
  # library(doParallel)
  # cl=makeCluster(18)
  # registerDoParallel(cl)
  
  library(samr)
  
  # if(missing(nbin)){nbin<-60}
  suppressPackageStartupMessages(library(SummarizedExperiment))
  suppressPackageStartupMessages(library(MultiAssayExperiment))
  
  dataset_name<-substr(dataset,1,nchar(dataset)-4)
  
  (dataset <- readRDS(paste0('data/',dataset)))
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
  if(dataset_name=='GSE62270-GPL17021'|dataset_name=='GSE78779-GPL17021'|dataset_name=='GSE81076-GPL18573'){
    dataset_gene_lstm=assays(dataset_gene)[['count']]
  }else{
    dataset_gene_lstm<-assays(dataset_gene)[['count']]
  }
  #dataset_gene_lstm<-assays(dataset_gene)[['count']]
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
  sname<-sname[sam==typeforFPR]
  #only contains the two selected cell group data
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  library(snow)
  library(pls)
  library(OGFSC)
  
  
  if(filter_method==1){
    filename<-paste0('Results/samr/',dataset_name,'-OGFSC-',nbin,'bin.csv')}else{
      if(filter_method==2){
        filename<-paste0('Results/samr/',dataset_name,'artifilt','.csv')
      }else{
        if(filter_method==3){
          filename<-paste0('Results/samr/',dataset_name,'without','.csv')
        }else{
          filename<-paste0('Results/samr/',dataset_name,'scmap','.csv')
        }
      }}
  
  
  snr_mean_m=NULL
  snr_cv_m=NULL
  snr_fz_m=NULL
  snr_var_m=NULL
  snr_mean=NULL
  snr_cv=NULL
  snr_fz=NULL
  snr_var=NULL
  
  #instances building & analysing
  ncells=dim(datarow)[2]
  maxsz<-floor(ncells*0.5)
  per=c(0.1,0.2,0.3,0.4)
  sort_gsz<-floor(ncells*per)
  
  #maxnum group 1 instance & analyse 
  
  set.seed(123)
  
  maxlis<-sample(1:ncells,2*maxsz,replace = F)
  
  instance_maxsz<-datarow[,maxlis]
  group_maxsz<-factor(c(rep('1',maxsz),rep('2',maxsz)))
  ##filter after sampling!!!
  if(filter_method==1){
    logmax<-log2(instance_maxsz+1)
    OGF = OGFSC(logmax,nBins = nbin, paral_option = 1, plot_option = 1,alpha = c(0.5,0.6,0.7,0.8,0.9,0.99,0.999))
    OGFSC_idx = OGF$OGFSC_idx
    instance_maxsz<-instance_maxsz[OGFSC_idx,] #gene selected after filtering by OGFSC
  }else{
    if(filter_method==2){
      if(dataset_name=='GSE81076-GPL18573'|dataset_name=='GSE62270-GPL17021'){
        instance_maxsz<-subset(instance_maxsz,rowSums(instance_maxsz>0)/dim(instance_maxsz)[2]>0.05)}else{
          instance_maxsz<-subset(instance_maxsz,rowSums(instance_maxsz>0)/dim(instance_maxsz)[2]>0.25)
        } ##article way to filter
    }else{
      if(filter_method==3){
        instance_maxsz<-subset(instance_maxsz,rowSums(instance_maxsz)>0)}else{
          source('C:/R project/DEmethods/scmap_filter.R')
          idx_ori=scmap_filter(datarow = instance_maxsz)
          instance_maxsz=instance_maxsz[idx_ori,]
        }
    }}
  #DE analyse
  # results <- apply(instance_maxsz, 1, ttestfun)
  ###results=p.adjust(results)   null data sets type I error control
  total_gene_num<-dim(instance_maxsz)[1]
  if(total_gene_num==0){
    FPR_ins_1<-NA
    cat('Warning: OGFSC filter with no genes left for max group!')
  }else{
    data=list(x=instance_maxsz,y=group_maxsz,geneid=as.character(1:nrow(instance_maxsz)),
              genenames=rownames(instance_maxsz),logged2=TRUE)
    res=try(samr(data,resp.type = 'Two class unpaired',nperms = 1000),silent = T)
    if('try-error'%in%class(res)){
      difgene_num_maxsz=NA
      FPR_ins_1=NA
      snr_mean_m=NA
      snr_cv_m=NA
      snr_fz_m=NA
      snr_var_m=NA
    }else{
    pv=samr.pvalues.from.perms(res$tt, res$ttstar)
    difgene_num_maxsz<-length(which(pv<0.05))
    FPR_ins_1<-difgene_num_maxsz/total_gene_num}
    results=pv
    idx1=which(p.adjust(pv,method = 'fdr')<=0.05)
    if(length(idx1)<=3){
      snr_mean_m=10086
      snr_cv_m=10086
      snr_fz_m=10086
      snr_var_m=10086
    }else{
      m_m=apply(instance_maxsz,1,logmean)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_mean_m=snr(m_s,m_ns)
      m_m=apply(instance_maxsz,1,logvar)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_var_m=snr(m_s,m_ns)
      m_m=apply(instance_maxsz,1,cv)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_cv_m=snr(m_s,m_ns)
      m_m=apply(instance_maxsz,1,fz)
      m_s=m_m[idx1]
      m_ns=m_m[-idx1]
      snr_fz_m=snr(m_s,m_ns)
    }
    cat('Complete maxmum group\n')}
  
  #make the other sizes' instance , replicate these instances 5 times, and do samr analyses
  
  repnum<-c(1:5)
  
  difgen_num<-array(NA,dim = c(1,5*length(per)))
  FPR_ge<-array(NA,dim = c(1,5*length(per)))
  n<--1
  for (j in sort_gsz) {
    ge_list<-c(rep(1,2*j))
    n<-n+1 
    for (i in repnum) {
      set.seed(i*j)
      ge_list<-sample(ncells,2*j)
      group_sz<-factor(c(rep('1',j),rep('2',j)))
      instance_ge<-datarow[,ge_list]####if sampling from instance_maxsz, OGFSC might have some problems
      if(filter_method==1){
        logge<-log2(instance_ge+1)
        OGF = OGFSC(logge,nBins = nbin, paral_option = 1, plot_option = 0,alpha = c(0.5,0.6,0.7,0.8,0.9,0.99,0.999))
        OGFSC_idx = OGF$OGFSC_idx
        instance_ge<-instance_ge[OGFSC_idx,] #gene selected after filtering by OGFSC
      }else{
        if(filter_method==2){
          if(dataset_name=='GSE81076-GPL18573'|dataset_name=='GSE62270-GPL17021'){
            instance_ge<-subset(instance_ge,rowSums(instance_ge>0)/dim(instance_ge)[2]>0.05)
          }else{
            instance_ge<-subset(instance_ge,rowSums(instance_ge>0)/dim(instance_ge)[2]>0.25)
          }
          ##article way to filter
        }else{
          if(filter_method==3){
            instance_ge<-subset(instance_ge,rowSums(instance_ge)>0)
          }else{
            source('C:/R project/DEmethods/scmap_filter.R')
            idx_ori=scmap_filter(datarow = instance_ge)
            instance_ge=instance_ge[idx_ori,]
          }
        }}
      
      ###results=p.adjust(results)
      total_gene_num<-dim(instance_ge)[1]
      if(total_gene_num==0){
        FPR_ge[1,n*5+i]<-NA
        cat('Warning: OGFSC filter with no genes left for group',j,',instance',i,'\n')
      }else{
        data=list(x=instance_ge,y=group_sz,geneid=as.character(1:nrow(instance_ge)),
                  genenames=rownames(instance_ge),logged2=TRUE)
        res=try(samr(data,resp.type = 'Two class unpaired',nperms = 1000),silent = T)
        if('try-error'%in%class(res)){
          difgen_num[1,n*5+i]=NA
          FPR_ge[1,n*5+i]=NA
          snr_mean[n*5+i]=NA
          snr_cv[n*5+i]=NA
          snr_fz[n*5+i]=NA
          snr_var[n*5+i]=NA
        }else{
          pv=samr.pvalues.from.perms(res$tt, res$ttstar)
          difgen_num[1,n*5+i]<- length(which(pv<0.05))
          FPR_ge[1,n*5+i]<-difgen_num[1,n*5+i]/total_gene_num 
          results=pv
          idx1=which(p.adjust(results,method = 'fdr')<=0.05)
          if(length(idx1)<=3){
            snr_mean[n*5+i]=10086
            snr_cv[n*5+i]=10086
            snr_fz[n*5+i]=10086
            snr_var[n*5+i]=10086
          }else{
            m_m=apply(instance_ge,1,logmean)
            m_s=m_m[idx1]
            m_ns=m_m[-idx1]
            snr_mean[n*5+i]=snr(m_s,m_ns)
            m_m=apply(instance_ge,1,logvar)
            m_s=m_m[idx1]
            m_ns=m_m[-idx1]
            snr_var[n*5+i]=snr(m_s,m_ns)
            m_m=apply(instance_ge,1,cv)
            m_s=m_m[idx1]
            m_ns=m_m[-idx1]
            snr_cv[n*5+i]=snr(m_s,m_ns)
            m_m=apply(instance_ge,1,fz)
            m_s=m_m[idx1]
            m_ns=m_m[-idx1]
            snr_fz[n*5+i]=snr(m_s,m_ns)
          }
        }
        
        cat('Complete groupsize of',j,',instacnce',i,'...Please wait\n',total_gene_num)}
    }
  }
  
  FPR<-c(FPR_ge,FPR_ins_1)
  group_size<-c(rep(sort_gsz,each=5),maxsz)
  differential_gene_num<-c(difgen_num,difgene_num_maxsz)
  SNR_cv=c(snr_cv,snr_cv_m)
  SNR_fz=c(snr_fz,snr_fz_m)
  SNR_mean=c(snr_mean,snr_mean_m)
  SNR_var=c(snr_var,snr_var_m)
  Results_dataframe<-data.frame(group_size,differential_gene_num,FPR,SNR_mean,SNR_var,SNR_cv,SNR_fz)
  #write the results to table
  
  write.csv(Results_dataframe,file = filename,quote = T)
  
  
  cat('Data analyse complete!')
  cat('\n ncells',ncells)
  
  return(Results_dataframe)
  ##stopCluster(cl)
}




