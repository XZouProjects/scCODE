#' Consensus Optimization of Differentially Expressed gene detection in single-cell RNA-seq data
#'
#' CODE presents personalized evaluation of different combination of gene filtering methods and DE methods for individual scRNA-seq data,
#' and provides with personalized optimal DE result.
#' @param datarow matrix, data matrix after filtering.
#' @param group vector, cell type information (eg: 1 for cell type 1, 2 for cell type 2)
#' @param light True or False, run CODE in a light version (apply on part of methods which averagely perform well), default as False (runing with all methods).
#' @return AUCC, matrix results of combinations of methods.
#' @return CDO, matrix results of combinations of methods.
#' @return Optimal solution, optimal combination of methods for the data.
#' @return Optimal results, optimal DE results detected by the optmal solution. ## DE gene name, p.adj, foldchange (cell type 1 vs celltype2, appeared firstly in group would be celltype1).
#' @author Jiawei Zou
#' @export
#' @examples
#' set.seed(123)
#' datarow=data_sample
#' group=group_sample
#' run_CODE(datarow,group,light=TRUE)
run_CODE<-function(datarow,group,light=TRUE){
  data_ori<-datarow
  group_ori<-group
  if(light==TRUE){
    DEmethods<-rev(c('MAST','t_test','BPSC','wilcox_test','DESeq2','limma','edgeR'))
    filtermethods<-c(1:3)
  }else{
    DEmethods<-rev(c('scDD','MAST','t_test','BPSC','wilcox_test','samr','DESeq2','limma','DEsingle','edgeR'))
    filtermethods<-c(1:4)
  }
  padj_list<-list()
  for(filt in filtermethods){
    if(filt==1){
      ###filtering by OGFSC
      idx<-CODE.filter_OGFSC(datarow = data_ori)
    }else{
      if(filt==2){
        ###filtering by conquer
        idx<-CODE.filter_conquer(data_ori)
      }else{
        if(filt==3){
          ###filtering by scmap
          idx<-CODE.filter_scmap(data_ori)
        }else{
          ###without filtering
          idx=which(rowSums(data_ori)>0)
        }
      }
    }

  datarow<-data_ori[idx,]
  idx_c=which(colSums(datarow)==0)
  if(length(idx_c)!=0){
    datarow=datarow[,-idx_c]
    group=group_ori[-idx_c]
  }
  ###gene name
  genename<-rownames(datarow)
  list_temp<-list()
  for (de in DEmethods) {
    run_func<-paste0('CODE.',de,'(datarow,group)')
    res_temp<-eval(parse(text = run_func))
    ### adjust p value<=0.05
    res_temp<-stats::p.adjust(res_temp,method = 'fdr')
    idx_sig=which(res_temp<=.05)
    ## top 500 DE genes
    n_gene=min(500,length(idx_sig))
    idx_sig=order(res_temp)[1:n_gene]
    ### padj record
    padj_temp<-res_temp[idx_sig]
    res_temp<-genename[idx_sig]
    ### DE gene record
    list_temp=c(list_temp,list(res_temp))
    padj_list=c(padj_list,list(padj_temp))

  }
  if(filt==1){
    list_OGFSC<-list_temp
  }else{
    if(filt==2){
      list_conquer<-list_temp
    }else{
      if(filt==3){
        list_scmap<-list_temp
      }else{
        list_without<-list_temp
      }
    }
  }
  }
  if(light==F){
    list_all<-c(list_OGFSC,list_conquer,list_scmap,list_without)
    ##CDO
    CDO_results<-cal_cdo(list_all)
    rescdo<-matrix(CDO_results,nrow = length(filtermethods),byrow = T)
    rownames(rescdo)=c('OGFSC','conquer','scmap','without_filter')
    colnames(rescdo)=DEmethods
    ##AUCC
    aucc1=cal_aucc(list_OGFSC)
    aucc2=cal_aucc(list_conquer)
    aucc3=cal_aucc(list_scmap)
    aucc4=cal_aucc(list_without)
    ##Average AUCC
    o1=(colSums(aucc1)-1)/(length(DEmethods)-1)
    o2=(colSums(aucc2)-1)/(length(DEmethods)-1)
    o3=(colSums(aucc3)-1)/(length(DEmethods)-1)
    o4=(colSums(aucc4)-1)/(length(DEmethods)-1)
    resaucc=data.frame(rbind(o1,o2,o3,o4))
    rownames(resaucc)=c('OGFSC','conquer','scmap','without_filter')
    colnames(resaucc)=DEmethods
  }else{
    list_all<-c(list_OGFSC,list_conquer,list_scmap)
    ##CDO
    CDO_results<-cal_cdo(list_all)
    rescdo<-matrix(CDO_results,nrow = length(filtermethods),byrow = T)
    rownames(rescdo)=c('OGFSC','conquer','scmap')
    colnames(rescdo)=DEmethods
    ##AUCC
    aucc1=cal_aucc(list_OGFSC)
    aucc2=cal_aucc(list_conquer)
    aucc3=cal_aucc(list_scmap)
    ##Average AUCC
    o1=(colSums(aucc1)-1)/(length(DEmethods)-1)
    o2=(colSums(aucc2)-1)/(length(DEmethods)-1)
    o3=(colSums(aucc3)-1)/(length(DEmethods)-1)
    resaucc=data.frame(rbind(o1,o2,o3))
    rownames(resaucc)=c('OGFSC','conquer','scmap')
    colnames(resaucc)=DEmethods
  }
  res_all=rescdo+resaucc
  res_t=NULL
  for (nc in 1:nrow(res_all)) {
    res_t<-c(res_t,as.numeric(res_all[nc,]))
  }
  idx_optimal=which(res_t==max(res_t))
  ##optimal filter method
  optimal_filter<-c('OGFSC','conquer','scmap','without_filter')[idx_optimal%/%length(DEmethods)+1]
  ##optimal DE method
  if(idx_optimal%%length(DEmethods)==0){
    optimal_de<-DEmethods[length(DEmethods)]
  }else{
    optimal_de<-DEmethods[idx_optimal%%length(DEmethods)]
  }
  ##optimal solution exhibition
  optimal_solution<-paste0('The optimal solution for your data is ',optimal_filter,' + ',optimal_de)
  ## p.adj for optimal DE genes
  padj_optimal=padj_list[[idx_optimal]]
  DE_genes_opotimal=list_all[[idx_optimal]]
  idx_m=match(DE_genes_opotimal,rownames(data_ori))
  data_m=data_ori[idx_m,]
  ### appeared later in group as cell type 2
  idx2<-which(group_ori==unique(group_ori)[2])
  ##cal_fc for each optimal DE genes detected
  cal_fc<-function(data){
    ex_1<-data[-idx2]
    ex_2<-data[idx2]
    fc<-mean(ex_1)/mean(ex_2)
    return(fc)
  }
  fc_optimal<-apply(data_m, 1, cal_fc)
  Optimal_results<-cbind(DE_genes_opotimal,padj_optimal,fc_optimal)
  colnames(Optimal_results)=c('DE genes','p-adjust','fold-change')
  Ressults_all<-list()
  Ressults_all<-c(Ressults_all,list(resaucc),list(rescdo),list(optimal_solution),list(Optimal_results))
  names(Ressults_all)<-c('AUCC','CDO','Optimal_solution','Optimal_results')
  cat('\nAnalysis Finished\n')
  return(Ressults_all)
}

#' Calculate AUCC value for data
#' @param lista_same_filter list, top 500 DE genes detected by every DE method, calculated within the same filtering method group.
#' @return AUCC (Area Under Consistency Curve) results.
cal_aucc=function(lista_same_filter){
  leh=length(lista_same_filter)
  nam=names(lista_same_filter)
  res=matrix(nrow = leh,ncol = leh)
  for (i in 1:(leh-1)) {
    res[i,i]=1
    for (j in (i+1):leh) {
      l=min(length(lista_same_filter[[i]]),length(lista_same_filter[[j]]))
      a=lista_same_filter[[i]][1:l]
      b=lista_same_filter[[j]][1:l]
      aucc_n=0
      t0=0
      t1=0
      for (n in 1:l) {
        t1=2*n-length(unique(c(a[1:n],b[1:n])))
        tn=t1-t0
        if(tn==0){
          aucc_n=aucc_n+0
        }else{
          aucc_n=aucc_n+1/2+(l-n)*tn
        }
        t0=t1
      }
      res[i,j]=aucc_n/(l*l/2)
      #res[i,j]=length(intersect(lista_same_filter[[i]],lista_same_filter[[j]]))/length(union(lista_same_filter[[i]],lista_same_filter[[j]]))###intersect/union
      res[j,i]=res[i,j]
    }
  }
  res[leh,leh]=1
  rownames(res)=nam
  colnames(res)=nam
  return(res)
}

#' Calculate CDO values for data
#' @param listall list, top 500 DE genes detected by every combination of methods.
#' @return CDO (Consistent DE gene Order) result for every method
cal_cdo=function(listall){
  l=length(listall)
  temp=listall[[1]]
  for (i in 2:l) {
    temp=c(temp,listall[[i]])
  }
  test=table(temp)
  test=as.data.frame(test)
  test=test[order(test$Freq,decreasing = T),]
  ## detected by more than 90% of methods (can adjust to lower), and no less than 5 genes
  thres=0.9
  res=test$temp[which(test$Freq>=(l*thres))]
  while(length(res)<5){
    thres=thres-0.05
    res=test$temp[which(test$Freq>=(l*thres))]
  }
  pos_all=NULL
  if(length(res)<5){
    pos_all=rep(NA,l)
  }else{
    for (i in 1:l) {
      # temp_m=length(res)
      # temp_n=length(listall[[i]])
      # temp_max=(1+temp_m)*temp_m/2
      # temp_min=temp_m*(temp_n+1)
      idx_tm=match(res,listall[[i]])
      ###balance to be 0-1
      idx_tm[is.na(idx_tm)]=length(listall[[i]])+length(res)/2
      idx_ave=mean(idx_tm)-length(res)/2
      porp_ave=idx_ave/length(listall[[i]])
      porp_ave=1-porp_ave#####近似最大假设为1，此时认为prop_ave值越大表示detective的能力越强,近似归一化
      pos_all=c(pos_all,porp_ave)

    }
  }

  return(pos_all)

}
