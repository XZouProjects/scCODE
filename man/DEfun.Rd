% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_functions.R
\name{DEfun}
\alias{DEfun}
\alias{scCODE.BPSC}
\alias{scCODE.t_test}
\alias{scCODE.DESeq2}
\alias{scCODE.MAST}
\alias{scCODE.wilcox_test}
\alias{scCODE.DEsingle}
\alias{scCODE.edgeR}
\alias{scCODE.samr}
\alias{scCODE.limma}
\alias{scCODE.scDD}
\title{DE methods used in CODE}
\usage{
scCODE.BPSC(datarow, group)

scCODE.t_test(datarow, group)

scCODE.DESeq2(datarow, group)

scCODE.MAST(datarow, group)

scCODE.wilcox_test(datarow, group)

scCODE.DEsingle(datarow, group)

scCODE.edgeR(datarow, group)

scCODE.samr(datarow, group)

scCODE.limma(datarow, group)

scCODE.scDD(datarow, group)
}
\arguments{
\item{datarow}{matrix, data matrix after filtering.}

\item{group}{factor, cell group information for data}
}
\value{
results, the p value for each gene.
}
\description{
DE methods used to detect DE genes for scRNA-seq data in CODE (Consensus Optimization of Differentially Expressed gene detection).
}
\examples{
set.seed(123)
datarow=data_sample
group=group_sample
scCODE.t_test(datarow,group)
}
\author{
Jiawei Zou
}
