
realdata_snr_anlys<-function(typ,dataset,celltypenames,filter_method,nbin){
  
  
  if(missing(nbin)){nbin<-50}
  library(doParallel)
  # cl=makeCluster(8)
  # registerDoParallel(cl)
  
  
  # if(missing(nbin)){nbin<-60}
  suppressPackageStartupMessages(library(SummarizedExperiment))
  suppressPackageStartupMessages(library(MultiAssayExperiment))
  
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
    typeforFPR2=levels(cop)[4]
  }else{
    if(dataset_name=='GSE60749-GPL13112'|dataset_name=='GSE100911'|dataset_name=='GSE94383'){
      cop<-anota$characteristics_ch1.1
      sam=paste0(sam,cop)
    }else{
      cop<-anota$characteristics_ch1
    }
  }
  if(dataset_name=='EMTAB3929'){
    cop=anota$Characteristics.developmental.stage.
    sam=paste0(cop,sam)
  }
  if(dataset_name=='GSE62270-GPL17021'|dataset_name=='GSE78779-GPL17021'|dataset_name=='GSE81076-GPL18573'){
    dataset_gene_lstm=assays(dataset_gene)[['count']]
  }else{
    dataset_gene_lstm<-assays(dataset_gene)[['count']]
  }
  if(dataset_name=='GSE77847'){
    sam=cop
  }
  if(dataset_name=="GSE81076-GPL18573"){
    typeforFPR='cell type: CD13+ sorted cells'
    typeforFPR2="cell type: TGFBR3+ sorted cells"
  }
  #dataset_gene_lstm<-assays(dataset_gene)[['TPM']]
  data_matrix<-dataset_gene_lstm[,sam==celltypenames[1]|sam==celltypenames[2]]
  cop<-cop[sam==celltypenames[1]|sam==celltypenames[2]]
  sam<-sam[sam==celltypenames[1]|sam==celltypenames[2]]
  sname<-colnames(data_matrix)
  #all zero gene kick out
  data_matrix<-subset(data_matrix,rowSums(data_matrix)>0)
  if(typ==3){
    colnames(data_matrix)<-cop
    datarow<-data_matrix[,cop==typeforFPR|cop==typeforFPR2]
    cop=cop[cop==typeforFPR|cop==typeforFPR2]
    sam=cop
  }else{
    
    colnames(data_matrix)<-sam
    datarow<-data_matrix# data for FPR
  }
  # sname<-sname[sam==typeforFPR]
  #only contains the two selected cell group data
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  library(snow)
  library(pls)
  library(OGFSC)
  
  if(filter_method==1){
    filename<-paste0('./',dataset_name,'/OGFSC1/',nbin,'-bin')}else{
      if(filter_method==2){
        filename<-paste0('./',dataset_name,'/conquer1/')
      }
      if(filter_method==3){filename<-paste0('./',dataset_name,'/nofilter1/')}
      if(filter_method==4){filename<-paste0('./',dataset_name,'/scmap1/')}
    }
  
  
  #instances building & analysing
  # snr_mean=NULL
  # snr_cv=NULL
  # snr_fz=NULL
  # snr_var=NULL
  data_ori=datarow
  
  for (r in seq(10)) {
    id_de=0
    
    idx2=which(sam==unique(sam)[2])
    group_maxsz<-c(rep('1',dim(datarow)[2]))
    group_maxsz[idx2]='2'
    group_maxsz=factor(group_maxsz)
    idx1=which(sam!=unique(sam)[2])
    set.seed(r+2021)
    r1=sample(idx1,length(idx1),replace = T)
    r2=sample(idx2,length(idx2),replace = T)
    
    datarow=data_ori[,c(r1,r2)]
    group_maxsz=group_maxsz[c(r1,r2)]
  
    snr_mean=NULL
    snr_cv=NULL
    snr_fz=NULL
    snr_var=NULL
  
  
  
  ####remove null cell(always for UMI data)
  
  tot0=dim(datarow)[1]
  
  ##filter after sampling!!!
  if(filter_method==1){
    logmax<-log2(datarow+1)
    OGF = OGFSC(logmax,nBins = nbin, paral_option = 1, plot_option = 1,alpha = c(0.5,0.6,0.7,0.8,0.9,0.99,0.999))
    idx_ori = OGF$OGFSC_idx
    #gene selected after filtering by OGFSC
  }else{
    if(filter_method==2){
      if(dataset_name=='GSE81076-GPL18573'|dataset_name=='GSE62270-GPL17021'){
        idx_ori<-which(rowSums(datarow>0)/dim(datarow)[2]>0.05)}else{
          idx_ori=which(rowSums(datarow>0)/dim(datarow)[2]>0.25)
        } 
      ##conquer
    }else{
      if(filter_method==3){idx_ori=which(rowSums(datarow)>0)}
      if(filter_method==4){
        source('C:/R project/DEmethods/scmap_filter.R')
        idx_ori=scmap_filter(datarow = datarow)}
    }
  }
  #DE analyse
  filter_idx=idx_ori
  fil_f=rep(0,tot0)
  fil_f[filter_idx]=1### =1 represent unfiltered genes
  
  mean_tgroup=function(expression_simulate){
    l=length(expression_simulate)
    idx_l=idx2
    a=mean(expression_simulate[-idx_l])
    b=mean(expression_simulate[idx_l])
    return(c(a,b))
  }
  expression_group=t(apply(datarow, 1, mean_tgroup))
  
  times_detected=rep(0,tot0)
  
  data_gene_informaion=cbind(expression_group,fil_f,times_detected)
  
  colnames(data_gene_informaion)=c('Mean_A','Mean_B','in_Filtered','Times_detected')
  
  datarow=datarow[idx_ori,]
  #DE analyse
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
  #######
  ##111111
  #ttest
  
  ttestfun<-function(x){
    tx<-x[idx2]
    ty<-x[-idx2]
    result<-t.test(tx,ty)
    pvalue<-result$p.value
    return(pvalue)
  }
  
  
  results <- apply(logdata, 1, ttestfun)
  ###results=p.adjust(results)   null data sets type I error control
  if(total_gene_num==0){
    FPR_ins_1<-NA
    cat('Warning: OGFSC filter with no genes left for max group!')
  }else{
    denum=min(length(which(results<=.05)))
    idx=order(results)[1:denum]
    denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
    DE_p<-length(idx)/total_gene_num
    DE_p1<-length(idx1)/total_gene_num
    DEgene=rnames[idx]
    DEgene1=rnames[idx1]
    
    regulate1=c(0,0,regulate[idx])
    regulate2=c(0,0,regulate[idx1])
    de_results=c(ncells,DE_p,DEgene)
    de_results1=c(ncells,DE_p1,DEgene1)
    
    DE_results=data.frame(DEgene=de_results,regu=regulate1)
    DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
    
    write.csv(DE_results,paste0(filename,'-ttest-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-ttest-adj_results-round',r,'.csv'))
    ### overlap
    denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
    idx_n=idx_ori[idx1]
    data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
    ###
    ##
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
  #######
  #222222
  ####BPSC
  library(BPSC)
  cl=makeCluster(18)
  registerDoParallel(cl)
  
  
  controlIds=which(group_maxsz==1)
  design=model.matrix(~group_maxsz)
  coef=2
  res=BPglm(data=datarow, controlIds=controlIds, design=design, coef=coef, estIntPar=F, useParallel=T)
  ###results=p.adjust(results)   null data sets type I error control
  total_gene_num<-dim(datarow)[1]
  if(total_gene_num==0){
    FPR_ins_1<-NA
    cat('Warning: OGFSC filter with no genes left for max group!')
  }else{
    results=res$PVAL
    denum=min(length(which(results<=.05)))
    idx=order(results)[1:denum]
    denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
    DE_p<-length(idx)/total_gene_num
    DE_p1<-length(idx1)/total_gene_num
    DEgene=rnames[idx]
    DEgene1=rnames[idx1]
    
    regulate1=c(0,0,regulate[idx])
    regulate2=c(0,0,regulate[idx1])
    de_results=c(ncells,DE_p,DEgene)
    de_results1=c(ncells,DE_p1,DEgene1)
    
    DE_results=data.frame(DEgene=de_results,regu=regulate1)
    DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
    write.csv(DE_results,paste0(filename,'-BPSC-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-BPSC-adj_results-round',r,'.csv'))
    stopCluster(cl)
    denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
    idx_n=idx_ori[idx1]
    data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
    ###
    ##
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
  
  ###
  ##33333
  ##DESeq2
  
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
  total_gene_num<-dim(datarow_in)[1]
  if(total_gene_num==0){
    cat('Warning: OGFSC filter with no genes left for max group!')
  }else{
    results=res_ins_maxsz$pvalue
    results=p.adjust(results,method = 'fdr')
    if(length(results)>length(total_gene_num)){
      results=results[-length(results)]
    }
    denum=min(length(which(results<=.05)))
    idx=order(results)[1:denum]
    denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
    DE_p<-length(idx)/total_gene_num
    DE_p1<-length(idx1)/total_gene_num
    DEgene=rnames[idx]
    DEgene1=rnames[idx1]
    
    regulate1=c(0,0,regulate[idx])
    regulate2=c(0,0,regulate[idx1])
    de_results=c(ncells,DE_p,DEgene)
    de_results1=c(ncells,DE_p1,DEgene1)
    
    DE_results=data.frame(DEgene=de_results,regu=regulate1)
    DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
    write.csv(DE_results,paste0(filename,'-DESeq2-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-DESeq2-adj_results-round',r,'.csv'))
    denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
    idx_n=idx_ori[idx1]
    data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
    ###
    ##
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
  
  ###MAST
  ##44444
  #
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
  total_gene_num<-dim(datarow)[1]
  if(total_gene_num==0){
    FPR_ins_1<-NA
    cat('Warning: OGFSC filter with no genes left for max group!')
  }else{
    results=fcHurdle$`Pr(>Chisq)`
    denum=min(length(which(results<=.05)))
    idx=order(results)[1:denum]
    denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
    DE_p<-length(idx)/total_gene_num
    DE_p1<-length(idx1)/total_gene_num
    DEgene=rnames[idx]
    DEgene1=rnames[idx1]
    
    regulate1=c(0,0,regulate[idx])
    regulate2=c(0,0,regulate[idx1])
    de_results=c(ncells,DE_p,DEgene)
    de_results1=c(ncells,DE_p1,DEgene1)
    
    DE_results=data.frame(DEgene=de_results,regu=regulate1)
    DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
    write.csv(DE_results,paste0(filename,'-MAST-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-MAST-adj_results-round',r,'.csv'))
    denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
    idx_n=idx_ori[idx1]
    data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
    ###
    ##
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
  
  
  ###wilcox test
  ###555555
  ##
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
  if(total_gene_num==0){
    FPR_ins_1<-NA
    cat('Warning: OGFSC filter with no genes left for max group!')
  }else{
    
    denum=min(length(which(results<=.05)))
    idx=order(results)[1:denum]
    denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
    DE_p<-length(idx)/total_gene_num
    DE_p1<-length(idx1)/total_gene_num
    DEgene=rnames[idx]
    DEgene1=rnames[idx1]
    
    regulate1=c(0,0,regulate[idx])
    regulate2=c(0,0,regulate[idx1])
    de_results=c(ncells,DE_p,DEgene)
    de_results1=c(ncells,DE_p1,DEgene1)
    
    DE_results=data.frame(DEgene=de_results,regu=regulate1)
    DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
    write.csv(DE_results,paste0(filename,'-wilcoxtest-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-wilcoxtest-adj_results-round',r,'.csv'))
    denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
    idx_n=idx_ori[idx1]
    data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
    ###
    ##
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
  
  ###DEsingle
  ###66666
  #
  library(DEsingle)
  results <- try(DEsingle(counts = datarow,group = group_maxsz,parallel = T),silent = T)
  if('try-error'%in% class(results)){
    DE_results=NULL
    DE_results1=NULL
    write.csv(DE_results,paste0(filename,'-DEsingle-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-DEsingle-adj_results-round',r,'.csv'))
  }else{
    results=results$pvalue
    
    denum=min(length(which(results<=.05)))
    idx=order(results)[1:denum]
    denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
    DE_p<-length(idx)/total_gene_num
    DE_p1<-length(idx1)/total_gene_num
    DEgene=rnames[idx]
    DEgene1=rnames[idx1]
    
    regulate1=c(0,0,regulate[idx])
    regulate2=c(0,0,regulate[idx1])
    de_results=c(ncells,DE_p,DEgene)
    de_results1=c(ncells,DE_p1,DEgene1)
    
    DE_results=data.frame(DEgene=de_results,regu=regulate1)
    DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
    write.csv(DE_results,paste0(filename,'-DEsingle-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-DEsingle-adj_results-round',r,'.csv'))
    denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
    idx_n=idx_ori[idx1]
    data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
    ###
    ##
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
  
  # ##edgeR
  # #7777777
  # #
  cl=makeCluster(16)
  registerDoParallel(cl)
  library(edgeR)
  
  dgelist =DGEList(counts = datarow, group = group_maxsz)
  dgelist_norm <- try(calcNormFactors(dgelist, method = 'TMM'),silent=T)
  if('try-error'%in%class(dgelist_norm)){
    DE_results=NULL
    DE_results1=NULL
    write.csv(DE_results,paste0(filename,'-edegeR-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-edgeR-adj_results-round',r,'.csv'))
  }else{
    design=model.matrix(~group_maxsz)
    dge <- estimateDisp(dgelist_norm, design, robust = F) #dispersion not robust
    fit <- glmFit(dge, design, robust = F)     #fit
    lrt <- glmLRT(fit)   #test
    results=lrt$table$PValue
    
    denum=min(length(which(results<=.05)))
    idx=order(results)[1:denum]
    denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
    DE_p<-length(idx)/total_gene_num
    DE_p1<-length(idx1)/total_gene_num
    DEgene=rnames[idx]
    DEgene1=rnames[idx1]
    
    regulate1=c(0,0,regulate[idx])
    regulate2=c(0,0,regulate[idx1])
    de_results=c(ncells,DE_p,DEgene)
    de_results1=c(ncells,DE_p1,DEgene1)
    
    DE_results=data.frame(DEgene=de_results,regu=regulate1)
    DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
    write.csv(DE_results,paste0(filename,'-edgeR-results-round',r,'.csv'))
    write.csv(DE_results1,paste0(filename,'-edgeR-adj_results-round',r,'.csv'))
    denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
    idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
    idx_n=idx_ori[idx1]
    data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
  }
  stopCluster(cl)
  ###
  ##
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
  
  # ###999999
  # ##
  # #samr
  library(samr)
  data=list(x=datarow,y=group_maxsz,geneid=as.character(1:nrow(datarow)),
            genenames=rownames(datarow),logged2=TRUE)
  res=samr(data,resp.type = 'Two class unpaired',nperms = 1000)
  results=samr.pvalues.from.perms(res$tt, res$ttstar)
  denum=min(length(which(results<=.05)))
  idx=order(results)[1:denum]
  denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
  idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
  DE_p<-length(idx)/total_gene_num
  DE_p1<-length(idx1)/total_gene_num
  DEgene=rnames[idx]
  DEgene1=rnames[idx1]
  
  regulate1=c(0,0,regulate[idx])
  regulate2=c(0,0,regulate[idx1])
  de_results=c(ncells,DE_p,DEgene)
  de_results1=c(ncells,DE_p1,DEgene1)
  
  DE_results=data.frame(DEgene=de_results,regu=regulate1)
  DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
  write.csv(DE_results,paste0(filename,'-samr-results-round',r,'.csv'))
  write.csv(DE_results1,paste0(filename,'-samr-adj_results-round',r,'.csv'))
  denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
  idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
  idx_n=idx_ori[idx1]
  data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
  ###
  ##
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
  ###
  ##10101010
  #limma
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
  
  denum=min(length(which(results<=.05)))
  idx=order(results)[1:denum]
  denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
  idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
  DE_p<-length(idx)/total_gene_num
  DE_p1<-length(idx1)/total_gene_num
  DEgene=rnames[idx]
  DEgene1=rnames[idx1]
  
  regulate1=c(0,0,regulate[idx])
  regulate2=c(0,0,regulate[idx1])
  de_results=c(ncells,DE_p,DEgene)
  de_results1=c(ncells,DE_p1,DEgene1)
  
  DE_results=data.frame(DEgene=de_results,regu=regulate1)
  DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
  write.csv(DE_results,paste0(filename,'-limma-results-round',r,'.csv'))
  write.csv(DE_results1,paste0(filename,'-limma-adj_results-round',r,'.csv'))
  denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
  idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
  idx_n=idx_ori[idx1]
  data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
  ###
  ##
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
  ###scDD 11
  library(scDD)
  library(SingleCellExperiment)
  names(group_maxsz)=colnames(datarow)
  condition=group_maxsz
  
  
  sce=SingleCellExperiment(assays=list(counts=datarow),colData=data.frame(condition))
  
  scDatEx.scran <- preprocess(sce, zero.thresh=1.0, scran_norm=TRUE)
  
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  scDatExSim <- scDD(scDatEx.scran, prior_param=prior_param, categorize = F,testZeroes=FALSE)
  
  res_scDD=results(scDatExSim)
  res_scDD=res_scDD$nonzero.pvalue
  results=res_scDD
  denum=min(length(which(results<=.05)))
  idx=order(results)[1:denum]
  denum1=min(length(which(p.adjust(results,method = 'fdr')<=.05)))
  idx1=order(p.adjust(results,method = 'fdr'))[1:denum1]
  DE_p<-length(idx)/total_gene_num
  DE_p1<-length(idx1)/total_gene_num
  DEgene=rnames[idx]
  DEgene1=rnames[idx1]
  
  regulate1=c(0,0,regulate[idx])
  regulate2=c(0,0,regulate[idx1])
  de_results=c(ncells,DE_p,DEgene)
  de_results1=c(ncells,DE_p1,DEgene1)
  
  DE_results=data.frame(DEgene=de_results,regu=regulate1)
  DE_results1=data.frame(DEgene=de_results1,regu=regulate2)
  write.csv(DE_results,paste0(filename,'-scDD-results-round',r,'.csv'))
  write.csv(DE_results1,paste0(filename,'-scDD-adj_results-round',r,'.csv'))
  denum_ori=min(length(which(p.adjust(results,method = 'fdr')<=.05)),500)
  idx1=order(p.adjust(results,method = 'fdr'))[1:denum_ori]
  idx_n=idx_ori[idx1]
  data_gene_informaion[idx_n,4]=data_gene_informaion[idx_n,4]+1
  cat('complete analysing dataset x!\n')
  cat('\n ncells',ncells)
  ###
  ##
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
  
  SNR_res=cbind(snr_var,snr_cv,snr_mean,snr_fz)
  rownames(SNR_res)=c('ttest','BPSC','DESeq2','MAST','wilcoxtest','DEsingle','edgeR','samr','limma','scDD')
  write.csv(SNR_res,paste0(filename,'-SNR-round',r,'.csv'))
  write.csv(data_gene_informaion,paste0(filename,'-DEinfor-round',r,'.csv'))
  ##stopCluster(cl)
}
}



