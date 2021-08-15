#' Gene filtering methods used in CODE
#'
#' Gene filtering methods used in CODE (Consensus Optimization of Differentially Expressed gene detection).
#' @param datarow matrix, data matrix after filtering.
#' @return idx, the idx of remained genes.
#' @author Jiawei Zou
#' @examples
#' set.seed(123)
#' datarow=data_sample
#' CODE.filter_OGFSC(datarow)
#' @name filter
NULL

#' @rdname filter
#' @export
CODE.filter_OGFSC<-function(datarow){
  logmax<-log2(datarow+1)
  OGF = OGFSC::OGFSC(logmax,nBins = 50, paral_option = 0, plot_option = 1,alpha=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999))
  idx = OGF$OGFSC_idx #Genes remained after filtering by OGFSC
  return(idx)
}

#' @rdname filter
#' @export
CODE.filter_conquer<-function(datarow){
  threshold_conquer=0.25
  idx<-which(rowSums(datarow>0)/dim(datarow)[2]>=threshold_conquer)###conquer 0.25 threshold, -0.05 when cut too much genes
  while ((length(idx)<=2000)&(threshold_conquer>=0.1)) {
    threshold_conquer=threshold_conquer-0.05
    idx<-which(rowSums(datarow>0)/dim(datarow)[2]>=threshold_conquer)
  }
  return(idx)
}

#' @rdname filter
#' @export
CODE.filter_scmap<-function(datarow){
  sce=SingleCellExperiment::SingleCellExperiment(assays=list(counts=datarow))
  logcounts(sce) <- log2(counts(sce) + 1)
  SummarizedExperiment::rowData(sce)$feature_symbol <- rownames(sce)
  sce=scmap::selectFeatures(sce,n_features = 5000)
  idx <- which(as.data.frame(SummarizedExperiment::rowData(sce))$scmap_features)
  return(idx)
}
