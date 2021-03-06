#' CODE plot function
#'
#' Plot the results of run_CODE
#' @param list results of run_CODE
#' @return Heatmap of AUCC
#' @return Heatmap of CDO
#' @return Heatmap of Metrics (AUCC+CDO)
#' @author Jiawei Zou
#' @import pheatmap
#' @export
scCODE_plot<-function(list){
  ##AUCC
  AUCC<-as.matrix(list$AUCC)
  ##CDO
  CDO<-as.matrix(list$CDO)
  ##Metrics
  Metrics<-AUCC+CDO
  pheatmap::pheatmap(AUCC,main = 'Z-Score(AUCC)',cluster_rows = F,cluster_cols = F,legend_breaks = c(min(AUCC)+(max(AUCC)-min(AUCC))/4,max(AUCC)-(max(AUCC)-min(AUCC))/4),legend_labels = c('Unsuitable','Suitable'),fontsize = 13)
  pheatmap::pheatmap(CDO,main = 'Z-Score(CDO)',cluster_rows = F,cluster_cols = F,legend_breaks = c(min(CDO)+(max(CDO)-min(CDO))/4,max(CDO)-(max(CDO)-min(CDO))/4),legend_labels = c('Unsuitable','Suitable'),fontsize = 13)
  pheatmap::pheatmap(Metrics,main = 'Z-Score(AUCC)+Z-Score(CDO)',cluster_rows = F,cluster_cols = F,legend_breaks = c(min(Metrics)+(max(Metrics)-min(Metrics))/4,max(Metrics)-(max(Metrics)-min(Metrics))/4),legend_labels = c('Unsuitable','Suitable'),fontsize = 13)

}
