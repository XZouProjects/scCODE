#####simulate with fc sample from gamma shape 4, rate set as below
##
#

base=c(2,3,4,5,6,7,8,10)####E, fold change
rate=4/base##shape=4
dep=0.02
source('./auroc.R')
source('./DE_gamma.R')
simualtion=function(rate){
  for (p in rate) {
    DEfor_gamma(p,dep)
  }
}
simualtion(rate)

