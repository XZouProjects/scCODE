setwd('./')

##DE methods
de_methods=rev(c('scDD','MAST','ttest','BPSC','wilcoxtest','samr','DESeq2','limma','DEsingle','edgeR'))

##SNR functions
cv=function(a){
  res=var(a)/mean(a)
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

##data input for FPR evaluation
data_sets=read.csv('./datas.csv')
for (m in de_methods) {
for (l in c(4)) {
  # l=2
  if(l==1){
    
    filename=paste0('Results/',m,'/results.csv')
    
  }else {
    if(l==2){
      
      filename=paste0('Results/',m,'/results_nof.csv')
    }else{
      if(l==3){
        filename=paste0('Results/',m,'/results_ar.csv')
      }else{
        filename=paste0('Results/',m,'/results_scmap.csv')
      }
      
    }
  }
  filter_method=l
  results=NULL


##ttest
for (i in 1:nrow(data_sets)) {
  
#####OGFSC
  
celltypenames<-c(data_sets$celltype1[i],data_sets$celltype2[i])
dataset<-paste0(data_sets$Data[i],'.rds')
typ<-data_sets$typ[i]
typeforFPR<-data_sets$typeforfpr[i]###one cell type used for null datasets building
nbin<-data_sets$nbin[i]
source(paste0(m,'fun.R'))
test_r=paste0(m,'_anlys(typ,dataset,celltypenames,typeforFPR,filter_method,nbin)')

res=eval(parse(text = test_r))

results=rbind(results,res)
}

write.csv(results,filename)
}
}

