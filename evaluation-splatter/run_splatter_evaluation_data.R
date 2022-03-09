#####splatter simulate
require(splatter)
require(scater)
require(R.matlab)

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

setwd('/data')

source('./New_Metrics.R')
source('./DE_functions.R')
source('./Filter_functions.R')
source('./auroc.R')

# nts<-5
###simulated data selected

  dataset_name<-'EMTAB_h_k'
  
  data_new<-readRDS(paste0(dataset_name,'.rds'))
  datarow<-data_new[,which(colnames(data_new)==unique(colnames(data_new))[1])]
  #instances building & analysing
  ncells=dim(datarow)[2]
  maxsz<-floor(ncells*0.5)
  datarow<-datarow[which(rowSums(datarow)>0),]
  # dir.create(paste0('C:/R project/DEmethods/simulate buquan/',dataset_name,'/'))
  
  ######splatter
  
  params <- splatEstimate(datarow)
  
  params <- setParam(params, 'batchCells',2*ncells)
  
  
  
  for (j in c(0.05,0.1,0.15,0.2,0.25)) {
    
    #maxnum group 1 instance & analyse
    
    params <- setParam(params, 'de.prob',j)
    
    sim.groups <- splatSimulate(params,group.prob = c(0.5, 0.5), de.facLoc = 0.3, method = "groups",
                                verbose = FALSE)
    # sim.groups <- logNormCounts(sim.groups)
    geneinfo=rowData(sim.groups)
    des=rep(1,length(geneinfo$DEFacGroup1))
    
    des[which(geneinfo$DEFacGroup1==1&geneinfo$DEFacGroup2==1)]=0#####non_DE
    groupinfo=colData(sim.groups)$Group
    group=rep('1',length(groupinfo))
    group[which(groupinfo=='Group2')]='2'
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
    
    
    instance_maxsz<-counts(sim.groups)
    
    idx_non_zero<-which(rowSums(instance_maxsz)>0)
    
    instance_maxsz<-instance_maxsz[idx_non_zero,]
    
    des<-des[idx_non_zero]
    
    ####gene num in raw data
    tot0=dim(instance_maxsz)[1]
    
    DEAnum=sum(des)
    
    non_DEAnum<-tot0-DEAnum
    
    datatest=instance_maxsz
    
    genename_all=rownames(instance_maxsz)
    
    colnames(datatest)<-paste0('cell',seq(ncol(datatest)))
    
    
    
    data_ori<-datatest
    data_ori<-as.matrix(data_ori)
    # data_ori<-data_ori[which(rowSums(data_ori)>0),]
    group_ori<-group
    # nts<-stats::median(c(nts,5,10))### top 5 - 10 consensus
    
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
    ###gene record for CDO and AUCC
    list_record<-list()
    ###auroc for each method
    auroc_record<-NULL
    ###FC record for consensus selection
    list_fc<-list()
    ### DE gene record for Final results & jaacard
    list_ori<-list()
    ### p value adjust record for rank gene
    list_pvl=list()
    ### gene order for calculating CDO
    list_gene_order<-list()
    ###DE gene number detected by each methods record
    de_num_r=NULL
    
    #######auroc_record
    auroc_record<-NULL
    ########TPR record
    tpr_record<-NULL
    #######FDP record
    fdp_record<-NULL
    
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
    
    
    ###TPR-a,TPR,FPR for DE gene
    TPR=NULL
    TPR_A=NULL### AFTER filtering
    FPR_A=NULL
    FPR<-NULL
    ### For top N genes
    
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
      
      des1=des[idx]####des temp
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
        idx_aaa=order(res_temp)
        
        ########deso, reordered des, record for metrics calculating
        des_o=des1[idx_aaa]
        ### adjusted p value record
        list_pvl=c(list_pvl,list(res_temp))
        
        ## all gene order for each method, used to calculate CDO
        gene_rank=names(res_temp)[order(res_temp)]
        
        list_gene_order=c(list_gene_order,list(gene_rank))
        
        ##record for CDO and AUCC
        idx_sig=which(res_temp<=.05)
        
        ###record auroc
        auroc_o=auroc_c(des_o,1)
        auroc_record<-c(auroc_record,auroc_o)
        
        
        ########record gene for CDO and AUCC
        ##de num
        de_num_r=c(de_num_r,length(idx_sig))
        ## top DEAnum DE genes for calculating AUCC and CDO
        n_gene=min(DEAnum,length(idx_sig))
        idx_sig=order(res_temp)[1:n_gene]
        ### DEAnum DE gene detected record
        padj_temp<-res_temp[idx_sig]
        gene_DEAnum<-genename[idx_sig]
        ### Top DEAnum DE gene record
        list_record=c(list_record,list(gene_DEAnum))
        
        
        
        ###record all DE gene results
        idx_sig=which(res_temp<=.05)
        ###
        n_gene=length(idx_sig)
        
        #######tpr and fdp record
        tpr_record<-c(tpr_record,sum(des_o[1:n_gene])/sum(des_o))
        
        fdp_record<-c(fdp_record,1-sum(des_o[1:n_gene])/n_gene)
        
        idx_sig=order(res_temp)[1:n_gene]
        ###Padj record
        padj_temp<-res_temp[idx_sig]
        ### All DE gene
        gene_all<-genename[idx_sig]
        
        ### FC record for each result, idx_FC, match back to the origin data
        idx_fc<-match(gene_all,genename_all)
        
        fc_temp<-fc_all[idx_fc]
        
        list_fc<-c(list_fc,list(fc_temp))
        
        ###ALL DE gene record for final results
        list_ori<-c(list_ori,list(gene_all))
        
        ########SNR bias
        
        ######mean fz var cv fc
        
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
        
        
        ###record overall TPR FPR (based on all genes before filtering)
        TPR=c(TPR,sum(des_o[1:n_gene])/DEAnum)
        TPR_A=tpr_record
        FPR_A=c(FPR_A,sum(des_o[1:n_gene]==0)/sum(des_o==0))
        FPR=c(FPR,sum(des_o[1:n_gene]==0)/non_DEAnum)
        ###top DEAnum gene
        
      }
      
    }
   
    list_all<-list_record
    ######gene and FDP record seperated by filtering methods
    con_per<-NULL
    # for (flt in c(4:1)) {
    
    ##CDO seperately by filtering methods
    CDO_results<-scCODE_cdo(list_all,list_gene_order)
    # write.csv(CDO_results,paste0('C:/R project/DEmethods/simulate buquan/',dataset_name,'/CDO',j,'.csv'))
    rescdo<-matrix(CDO_results,nrow = length(filtermethods),byrow = T)
    rownames(rescdo)=c('OGFSC','conquer','scmap','without_filter')
    colnames(rescdo)=DEmethods
    ##AUCC
    ##AUCC
    aucc1=scCODE_aucc(list_all[1:10])
    aucc2=scCODE_aucc(list_all[11:20])
    aucc3=scCODE_aucc(list_all[21:30])
    aucc4=scCODE_aucc(list_all[31:40])
    
    o1=(colSums(aucc1)-1)/(length(DEmethods)-1)
    o2=(colSums(aucc2)-1)/(length(DEmethods)-1)
    o3=(colSums(aucc3)-1)/(length(DEmethods)-1)
    o4=(colSums(aucc4)-1)/(length(DEmethods)-1)
    
    AUCC_results<-c(o1,o2,o3,o4)
    
    resaucc=data.frame(rbind(o1,o2,o3,o4))
    rownames(resaucc)=c('OGFSC','conquer','scmap','without_filter')
    colnames(resaucc)=DEmethods
    
    res_all<-z_score(apply(rescdo,2,as.numeric))+z_score(apply(resaucc,2,as.numeric))
    res_t=NULL
    
    for (nc in 1:nrow(res_all)) {
      res_t<-c(res_t,as.numeric(res_all[nc,]))
    }
    for (nts in c(5,10,20,40)) {
      
      idx_optimal=order(res_t,decreasing = T)[1:nts]####top suitable strategies consensus selected by scCODE
      consensus_res<-NULL
      info_all<-NULL
      for (idx_o in idx_optimal) {
        gene_temp<-list_ori[[idx_o]]
        n_gene_temp<-length(gene_temp)
        info_all<-c(info_all,n_gene_temp)
        fc_temp<-list_fc[[idx_o]]
        consensus_temp<-cbind(gene_temp,fc_temp)
        consensus_res<-rbind(consensus_res,consensus_temp)
       
      consensus_res<-as.data.frame(consensus_res)
      ###gene freq
      counta<-table(consensus_res$gene_temp)
      counta<-sort(counta,decreasing = T)
      ##logFC
      idx_c<-match(names(counta),consensus_res$gene_temp)
      fca<-consensus_res$fc_temp[idx_c]
      logfc<-log2(as.numeric(fca))
     
      consensus_res<-data.frame(counta,logfc)
      if(dim(consensus_res)[2]==2){
        consensus_res<-cbind(rownames(consensus_res),consensus_res)
      }
      colnames(consensus_res)=c('Gene_name','Detected_times','logFC')
      ####rank by logFC
      idxf<-order(abs(consensus_res$logFC),decreasing = T)
      test<-consensus_res[idxf,]
      
      ###test conbined results
      test<-test[order(test$Detected_times,decreasing = T),]
  
      #####DE gene rank first by freq, then by fc, scCODE results
      scCODE_res<-test
      scCODE_res<-scCODE_res$Gene_name
      idx_sccode<-match(scCODE_res,genename_all)
      des_sccode<-des[idx_sccode]
      
      #####Bias of scCODE
      logmean_de<-logmean_all[idx_sccode]
      
      logmean_record<-c(mean(logmean_de),logmean_record)
      
      logmean_non<-logmean_all[-idx_sccode]
      
      snr_mean_temp=snr(logmean_de,logmean_non)
      
      snr_mean_record<-c(snr_mean_temp,snr_mean_record)
      
      logvar_de<-logvar_all[idx_sccode]
      
      logvar_record<-c(mean(logvar_de),logvar_record)
      
      logvar_non<-logvar_all[-idx_sccode]
      
      snr_var_temp=snr(logvar_de,logvar_non)
      
      snr_var_record<-c(snr_var_temp,snr_var_record)
      
      cv_de<-cv_all[idx_sccode]
      
      cv_record<-c(mean(cv_de),cv_record)
      
      cv_non<-cv_all[-idx_sccode]
      
      snr_cv_temp=snr(cv_de,cv_non)
      
      snr_cv_record<-c(snr_cv_temp,snr_cv_record)
      
      fz_de<-fz_all[idx_sccode]
      
      fz_record<-c(mean(fz_de),fz_record)
      
      fz_non<-fz_all[-idx_sccode]
      
      snr_fz_temp=snr(fz_de,fz_non)
      
      snr_fz_record<-c(snr_fz_temp,snr_fz_record)
      
      logfc_temp<-log2fc_all[idx_sccode]
      
      logfc_record<-c(mean(logfc_temp),logfc_record)
     
      common_gene<-NULL
      for (mo in idx_optimal) {
        genetemp=names(list_pvl[[mo]])
        common_gene=c(common_gene,genetemp)
      }
      common_gene<-unique(common_gene)
      oth_gene<-setdiff(common_gene,scCODE_res)
      
      ####the other_gene rank
      idx_oth_gene=match(oth_gene,genename_all)
      des_oth_gene=des[idx_oth_gene]
      ###overall auroc for scCODE
      des_sccode=c(des_sccode,des_oth_gene)
      ###record TPR FPR for scCODE
      ###record TPR FPR
      ############all DE gene summary
      n_median<-dim(test)[1]
      
      #######record scCODE results
      TPR=c(sum(des_sccode[1:n_median])/DEAnum,TPR)
      TPR_A=c(sum(des_sccode[1:n_median])/sum(des_sccode),TPR_A)
      FPR_A=c((sum(des_sccode[1:n_median]==0)/sum(des_sccode==0)),FPR_A)
      FPR=c((sum(des_sccode[1:n_median]==0)/DEAnum),FPR)
      ###top DEAnum gene
      
    }
   
    Evaluate=cbind(TPR,TPR_A,FPR_A,FPR,logmean_record,logvar_record,cv_record,fz_record,logfc_record,snr_mean_record,snr_var_record,snr_cv_record,snr_fz_record)
    rownames(Evaluate)=c('scCODE_40','scCODE_20','scCODE_10','scCODE_5',paste0('OGFSC-',DEmethods),paste0('conquer-',DEmethods),paste0('scmap-',DEmethods),paste0('without_filter-',DEmethods))
    
    dir.create(paste0('./',dataset_name))
    write.csv(Evaluate,paste0('./',dataset_name,'/Evaluate-scCODE-',j,'.csv'))
    
    metrics_on_simulate<-cbind(auroc_record,fdp_record,tpr_record,de_num_r,CDO_results,AUCC_results)
    rownames(metrics_on_simulate)<-c(paste0('OGFSC-',DEmethods),paste0('conquer-',DEmethods),paste0('scmap-',DEmethods),paste0('without_filter-',DEmethods))
    write.csv(metrics_on_simulate,paste0('./',dataset_name,'/Evaluate-metrics-',j,'.csv'))
    
    cat(paste0('/nAnalysis Finished/n dep-'),j)
    }
  }
  
