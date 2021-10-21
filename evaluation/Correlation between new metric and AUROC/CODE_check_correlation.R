base=c(6)####E, fold change
###2,4,3,6,9
rate=4/base##shape=4
dep=0.02
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
source('./auroc.R')
source('./CODE_correlation.R')####FOR AUCC CDO correlation
simualtion=function(rate){
  for (p in rate) {
    SNR_gamma(p,dep)####FOR SNR
  }
}
simualtion(rate)
