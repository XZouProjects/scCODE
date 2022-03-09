#' Calculate AUCC value for data
#' @param lista_same_filter list, top 500 DE genes detected by every DE method, calculated within the same filtering method group.
#' @return AUCC (Area Under Consistency Curve) results.
scCODE_aucc=function(lista_same_filter){
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
#' @param list_order, list, ALL gene rank for each analysis strategy.
#' @return CDO (Consistent DE gene Order) result for every method
scCODE_cdo=function(listall,list_order){
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
  while((length(res)<5)&(thres>=0.75)){
    thres=thres-0.05
    res=test$temp[which(test$Freq>=(l*thres))]
  }

  pos_all=NULL
  ###calculate CDO in all gene order
  for (i in 1:l) {
    idx_tm=match(res,list_order[[i]])

    if(length(idx_tm)==0){
      porp_ave=0.05
    }else{

    ###balance to be 0-1
    idx_tm[is.na(idx_tm)]=(2*length(list_order[[i]])-length(res)+1)/2
    idx_ave=sum(idx_tm)
    porp_ave=(idx_ave-(length(res)*(length(res)+1))/2)/(length(res)*(length(list_order[[i]])-length(res)))
    if((length(res)*(length(list_order[[i]])-length(res)))==0){
      porp_ave=0.5
    }
    porp_ave=1-porp_ave#####CDO
    
    
    if(porp_ave<=0){
      porp_ave=0.05###punishment for too less DE gene detected
    }
    }

    pos_all=c(pos_all,porp_ave)

  }
  return(pos_all)

}

#' Calculate Z score for data
#' @param data data, numberic verctor.
#' @return Zscore for data.
z_score=function(data){
  sd_data=stats::sd(data)
  average_data=mean(data)
  zs_data=(data-average_data)/sd_data
  return(zs_data)
}
