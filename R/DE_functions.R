#' DE methods used in CODE
#'
#' DE methods used to detect DE genes for scRNA-seq data in CODE (Consensus Optimization of Differentially Expressed gene detection).
#'
#' @param datarow matrix, data matrix after filtering.
#' @param group factor, cell group information for data
#' @return results, the p value for each gene.
#' @author Jiawei Zou
#' @examples
#' set.seed(123)
#' datarow=data_sample
#' group=group_sample
#' scCODE.t_test(datarow,group)
#' @name DEfun
NULL

#' @rdname DEfun
#' @export
scCODE.BPSC <- function(datarow, group) {
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_workers <- 2L
  } else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
  }
  cl=parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  controlIds=which(group==1)
  design=stats::model.matrix(~group)
  coef=2
  res=BPSC::BPglm(data=datarow, controlIds=controlIds, design=design, coef=coef, estIntPar=F, useParallel=T)
  results=res$PVAL
  parallel::stopCluster(cl)
  return(results)
}


#' @rdname DEfun
#' @export
scCODE.t_test<-function(datarow,group){
  idx2=which(group==unique(group)[2])
  ttestfun<-function(x){
    tx<-x[idx2]
    ty<-x[-idx2]
    result<-stats::t.test(tx,ty)
    pvalue<-result$p.value
    return(pvalue)
  }
  logdata=log2(datarow+1)
  results <- apply(logdata, 1, ttestfun)
  return(results)
}

#' @rdname DEfun
#' @export
scCODE.DESeq2<-function(datarow,group){
  ##integer
  datarow_in<-apply(datarow, 2, as.integer)

  dds_ins_maxsz<-DESeq2::DESeqDataSetFromMatrix(datarow_in,S4Vectors::DataFrame(group),~group)
  des_ins_maxsz<-try(DESeq2::DESeq(dds_ins_maxsz,parallel = F),silent = T)
  if('try-error'%in% class(des_ins_maxsz)){
    fakegene=round(colSums(datarow_in)/dim(datarow_in)[1])+1
    datarow_in=rbind(datarow_in,fakegene)
    dds_ins_maxsz<-DESeq2::DESeqDataSetFromMatrix(datarow_in,S4Vectors::DataFrame(group),~group)
    des_ins_maxsz<-try(DESeq2::DESeq(dds_ins_maxsz,parallel = F),silent = T)
  }
  res_ins_maxsz<-DESeq2::results(des_ins_maxsz)
  results=res_ins_maxsz$pvalue
  if(length(results)>nrow(datarow)){
    results=results[-length(results)]
  }
  return(results)
}

#' @rdname DEfun
#' @export
scCODE.MAST<-function(datarow,group){
  logdata=log2(datarow+1)
  cData<-data.frame(wellKey=paste0('C',1:dim(logdata)[2]))
  fData<-data.frame(primerid=rownames(logdata))
  colnames(logdata)=paste0('C',1:dim(logdata)[2])
  sca<-MAST::FromMatrix(logdata,cData,fData)
  SummarizedExperiment::colData(sca)$cond<-group
  zlmcond<-try(MAST::zlm(~cond,sca,parallel = T),silent = T)
  summarycond<-MAST::summary(zlmcond,doLRT='cond2')
  summaryDt <- summarycond$datatable
  fcHurdle <- summaryDt[summaryDt$contrast=='cond2' & summaryDt$component=='H',]#logFC coefficients
  results=fcHurdle$`Pr(>Chisq)`
  return(results)
}

#' @rdname DEfun
#' @export
scCODE.wilcox_test<-function(datarow,group){
  idx2=which(group==unique(group)[2])
  wilcoxtestfun<-function(x){
    tx<-x[idx2]
    ty<-x[-idx2]
    result<-stats::wilcox.test(tx,ty)
    pvalue<-result$p.value
    return(pvalue)
  }
  logdata=log2(datarow+1)
  results <- apply(logdata, 1, wilcoxtestfun)
  return(results)
}

#' @rdname DEfun
#' @export
scCODE.DEsingle<-function(datarow,group){
  genename=rownames(datarow)
  datatest=apply(datarow, 2, as.integer)
  rownames(datarow)=genename
  datarow=data.frame(datarow)
  results <- try(DEsingle::DEsingle(counts = datarow,group = group,parallel = T),silent = T)
  idxdes<-match(rownames(datarow),rownames(results))
  results=results[idxdes,]
  results=results$pvalue
  return(results)
}

#' @rdname DEfun
#' @export
scCODE.edgeR<-function(datarow,group){
  dgelist=edgeR::DGEList(counts = datarow, group = group)
  ### TMM method
  dgelist_norm <- try(edgeR::calcNormFactors(dgelist, method = 'TMM'),silent=T)
  design=stats::model.matrix(~group)
  dge <- edgeR::estimateDisp(dgelist_norm, design, robust = F) #dispersion not robust
  fit <- edgeR::glmFit(dge, design, robust = F)     #fit
  lrt <- edgeR::glmLRT(fit)   #test
  results=lrt$table$PValue
  return(results)
}

#' @rdname DEfun
#' @export
scCODE.samr<-function(datarow,group){
  data=list(x=datarow,y=group,geneid=as.character(1:nrow(datarow)),
            genenames=rownames(datarow),logged2=TRUE)
  res=try(samr::samr(data,resp.type = 'Two class unpaired',nperms = 1000),silent = T)
  results=samr::samr.pvalues.from.perms(res$tt, res$ttstar)
  return(results)
}

#' @rdname DEfun
#' @export
scCODE.limma<-function(datarow,group){
  idx2=which(group==unique(group)[2])
  limmagroup<-c(rep('control',dim(datarow)[2]))
  limmagroup[idx2]='fake'
  limmagroup=factor(limmagroup)
  names(limmagroup)=colnames(datarow)
  design=stats::model.matrix(~0+limmagroup)
  colnames(design)=levels(factor(limmagroup))
  rownames(design)=colnames(datarow)
  norm=limma::voom(datarow,design =design,plot = FALSE)
  fit=limma::lmFit(norm,design,method = 'ls')
  contrast=limma::makeContrasts('control-fake',levels = design)
  fit2=limma::contrasts.fit(fit,contrast)
  fit2=limma::eBayes(fit2)
  results=fit2$p.value
  return(results)
}

#' @rdname DEfun
#' @export
#' @import SingleCellExperiment
scCODE.scDD<-function(datarow,group){
  names(group)=colnames(datarow)
  condition=group
  idx=which(colSums(datarow)==0)
  if(length(idx)!=0){
    datarow=datarow[,-idx]
    condition=condition[-idx]
  }

  sce=SingleCellExperiment::SingleCellExperiment(assays=list(counts=datarow),colData=data.frame(condition))

  scDatEx.scran <- scDD::preprocess(sce, zero.thresh=1.0, scran_norm=TRUE)

  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  scDatExSim <- scDD::scDD(scDatEx.scran, prior_param=prior_param, categorize = F,testZeroes=FALSE)

  res_scDD=scDD::results(scDatExSim)
  results=res_scDD$nonzero.pvalue
  return(results)
}
