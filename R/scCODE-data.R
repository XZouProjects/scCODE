#' A scRNA-seq dataset containg activated and naive CD4+ T cells
#'
#' A scRNA-seq dataset containg activated and naive CD4+ T cells. 13045 genes x 462 cells.
#'
#' @format A conut matrix 13045 genes x 462 cells:
#' \describe{
#'   \item{row}{gene name}
#'   \item{column}{cell name}
#'   ...
#' }
#' @source \url{https://www.science.org/doi/full/10.1126/science.aah4115}
"data_sccode"

#' group information of dataset.
#'
#' A vector of 1 and 2, the group information of the cells. Here, 1 for activated CD4+ T cells, 2 for naive CD4+ T cells.
#'
#' @format A factor length 462
#' \describe{
#'   \item{factor 1}{activated CD4+ T cells}
#'   \item{factor 2}{naive CD4+ T cells}
#'   ...
#' }
#' @source \url{https://www.science.org/doi/full/10.1126/science.aah4115}
"group_sccode"
