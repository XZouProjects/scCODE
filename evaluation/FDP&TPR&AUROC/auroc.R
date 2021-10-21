auroc_c=function(label,positive_label){
  auroc_c=0
  for (i in 1:length(label)) {
    if(label[i]==positive_label){
      temp=length(which(label[1:i]!=positive_label))
      fpr=temp/length(which(label!=positive_label))
      auroc_c=auroc_c+(1-fpr)/(length(which(label==positive_label)))
    }
  }
  return(auroc_c)
}


