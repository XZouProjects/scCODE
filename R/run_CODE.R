#' Consensus Optimization of Differentially Expressed gene detection in single-cell RNA-seq data
#'
#' CODE presents personalized evaluation of different combination of gene filtering methods and DE methods for individual scRNA-seq data,
#' and provides with personalized optimal DE result.
#' @param data matrix, data matrix after filtering.
#' @param group vector, cell type information (eg: 1 for cell type 1, 2 for cell type 2)
#' @param light True or False, run CODE in a light version (apply on part of methods which averagely perform well), default as False (runing with all methods).
#' @param nts number of top-ranked strategies for consensus optimization, 5 =< nts <= 10.
#' @param outdir Path to save DE gene identification results.
#' @return AUCC, matrix results of combinations of methods.
#' @return CDO, matrix results of combinations of methods.
#' @return Optimal solution, optimal combination of methods for the data.
#' @return Optimal results, optimal DE results detected by the optmal solution. ## DE gene name, p.adj, foldchange (cell type 1 vs celltype2, appeared firstly in group would be celltype1).
#' @author Jiawei Zou
#' @export
#' @examples
#' set.seed(123)
#' run_CODE(data_sample,group_sample,light=TRUE,outdir='./',nts=5)
run_CODE<-function(data,group,light=TRUE,outdir,nts=5){
  data_ori<-data
  group_ori<-group
  nts<-stats::median(c(nts,5,10))### top 5 - 10 consensus
  ### appeared later in group as cell type 2
  idx2<-which(group==unique(group)[2])
  ##cal_fc for each optimal DE genes detected
  cal_fc<-function(data){
    ex_1<-data[-idx2]
    ex_2<-data[idx2]
    fc<-mean(ex_1)/mean(ex_2)
    return(fc)
  }
  ###FC was recorded as the value of origin data input
  fc_all<-apply(data_ori, 1, cal_fc)
  ###
  genename_all=rownames(data_ori)
  ###light version to save time
  if(light==TRUE){
    DEmethods<-rev(c('MAST','t_test','BPSC','wilcox_test','DESeq2','limma','edgeR'))
    filtermethods<-c(1:3)
  }else{
    DEmethods<-rev(c('scDD','MAST','t_test','BPSC','wilcox_test','samr','DESeq2','limma','DEsingle','edgeR'))
    filtermethods<-c(1:4)
  }

  ###top 500 gene record for finding seed genes in CDO and  calculating AUCC
  list_record<-list()
  ###FC record for consensus selection
  list_fc<-list()
  ### DE gene detected by each analysis strategy
  list_ori<-list()
  ### p value adjust record for rank gene
  list_pvl=list()
  ### gene order for calculating CDO
  list_gene_order<-list()

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
    ## all gene order for each method, used to calculate CDO
    gene_rank=names(res_temp)[order(res_temp)]

    list_gene_order=c(list_gene_order,list(gene_rank))

    ##record for CDO and AUCC
    idx_sig=which(res_temp<=.05)

    ## top 500 DE genes for calculating AUCC and CDO
    n_gene=min(500,length(idx_sig))
    idx_sig=order(res_temp)[1:n_gene]
    ### 500 DE gene detected record
    gene_500<-genename[idx_sig]
    ### Top 500 DE gene record
    list_record=c(list_record,list(gene_500))
    ###record all DE gene results
    idx_sig=which(res_temp<=.05)
    ###
    n_gene=length(idx_sig)
    idx_sig=order(res_temp)[1:n_gene]
    ###Padj record
    padj_temp<-res_temp[idx_sig]
    list_pvl=c(list_pvl,list(padj_temp))
    ### All DE gene
    gene_all<-genename[idx_sig]
    ### FC record for each result
    idx_fc<-match(gene_all,genename_all)
    fc_temp<-fc_all[idx_fc]
    list_fc<-c(list_fc,list(fc_temp))
    ###ALL DE gene record
    list_ori<-c(list_ori,list(gene_all))

  }
  }
  ###using top 500 to find seed gene and calculate AUCC
  list_all<-list_record
  ##CDO
  CDO_results<-scCODE_cdo(list_all,list_order = list_gene_order)
  ###Z score normalization
  CDO_results<-z_score(CDO_results)
  ###
  rescdo<-matrix(CDO_results,nrow = length(filtermethods),byrow = T)
  rownames(rescdo)=c('OGFSC','conquer','scmap','no filter')[filtermethods]
  colnames(rescdo)=DEmethods
  ##AUCC
  ##AUCC
  aucc1=scCODE_aucc(list_all[1:length(DEmethods)])
  aucc2=scCODE_aucc(list_all[(length(DEmethods)+1):(2*length(DEmethods))])
  aucc3=scCODE_aucc(list_all[(2*length(DEmethods)+1):(3*length(DEmethods))])
  o1=(colSums(aucc1)-1)/(length(DEmethods)-1)
  o2=(colSums(aucc2)-1)/(length(DEmethods)-1)
  o3=(colSums(aucc3)-1)/(length(DEmethods)-1)
  if(light==F){
    aucc4=scCODE_aucc(list_all[(3*length(DEmethods)+1):(4*length(DEmethods))])
    o4=(colSums(aucc4)-1)/(length(DEmethods)-1)
    ###Z-score
    AUCC_results=z_score(c(o1,o2,o3,o4))
  }else{
    ###Z-score
    AUCC_results=z_score(c(o1,o2,o3))
  }
  resaucc<-matrix(AUCC_results,nrow = length(filtermethods),byrow = T)
  rownames(resaucc)=c('OGFSC','conquer','scmap','no filter')[filtermethods]
  colnames(resaucc)=DEmethods

  res_all=rescdo+resaucc

  # Methods rank and consensus of sCODE results, consensus optimization of all the DE genes detected by selected analysis strategies
  ##
  ###
  res_t=NULL
  for (nc in 1:nrow(res_all)) {
    res_t<-c(res_t,as.numeric(res_all[nc,]))
  }
  idx_optimal=order(res_t,decreasing = T)[1:nts]####top suitable strategies consensus
  consensus_res<-NULL
  info_all<-NULL
  for (i in idx_optimal) {
    gene_temp<-list_ori[[i]]
    n_gene_temp<-length(gene_temp)
    fc_temp<-list_fc[[i]]
    padj_temp<-list_pvl[[i]]
    consensus_temp<-cbind(gene_temp,fc_temp,padj_temp)
    consensus_res<-rbind(consensus_res,consensus_temp)
    ##optimal filter method
    optimal_filter<-c('OGFSC','conquer','scmap','no filter')[ceiling(i/length(DEmethods))]
    ##optimal DE method
    if(i%%length(DEmethods)==0){
      optimal_de<-DEmethods[length(DEmethods)]
    }else{
      optimal_de<-DEmethods[i%%length(DEmethods)]
    }
    ##information of the method
    info_temp<-c(n_gene_temp,optimal_filter,optimal_de)
    info_all<-rbind.data.frame(info_all,info_temp)
  }
  colnames(info_all)<-c('n_DEgene','Filter','DE')
  n_median<-stats::median(info_all$n_DEgene)
  consensus_res<-as.data.frame(consensus_res)
  ###gene freq
  counta<-table(consensus_res$gene_temp)
  counta<-sort(counta,decreasing = T)
  ##logFC
  idx_c<-match(names(counta),consensus_res$gene_temp)
  fca<-consensus_res$fc_temp[idx_c]
  logfc<-log2(as.numeric(fca))
  ###P adjust (mininum)
  consensus_res<-consensus_res[order(consensus_res$padj_temp),]
  idx_p<-match(names(counta),consensus_res$gene_temp)
  p_adj<-consensus_res$padj_temp[idx_p]
  ##consensus results combined
  consensus_res<-data.frame(counta,logfc,p_adj)
  colnames(consensus_res)=c('Gene_name','Detected_times','logFC','P_adjust')
  ####rank by logFC
  idxf<-order(abs(consensus_res$logFC),decreasing = T)
  test<-consensus_res[idxf,]

  ###test conbined results
  test<-test[order(test$Detected_times,decreasing = T),]
  #####n gene rank first by freq, then by fc
  consensus_res<-test[1:n_median,]

  rownames(consensus_res)=seq(nrow(consensus_res))

  Ressults_all<-c(list(data.frame(resaucc)),list(data.frame(rescdo)),list(data.frame(info_all)),list(data.frame(consensus_res)))
  names(Ressults_all)<-c('AUCC','CDO','Strategy info','scCODE results')
  writexl::write_xlsx(Ressults_all, paste0(outdir,'/scCODE results.xlsx'))
  cat('\nAnalysis Finished\n')
  return(Ressults_all)
}
