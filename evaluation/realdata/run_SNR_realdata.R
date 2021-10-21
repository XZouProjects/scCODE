setwd('./')
source('./snr_e_realdata.R')
data_sets=read.csv('./qrdata.csv')
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
for (l in 1:4) {
  #l=4
  filter_method=l
  for (i in c(3,6,9)) {
    celltypenames<-c(data_sets$cell1[i],data_sets$cell2[i])
    dataset<-paste0(data_sets$Data[i],'.rds')
    typ<-data_sets$typ[i]
    nbin<-data_sets$nbins[i]
    realdata_snr_anlys(typ,dataset,celltypenames,filter_method,nbin)
  }
}

