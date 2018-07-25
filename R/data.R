#' Simulated bicluster datasets
#' 
#' Simulated matrices containing biclusters and mild normally-distributed noise.
#' Generated using \code{genSimData}.
#' 
#' @format A list of 4 matrices named:
#' \describe{
#'  \item{nonoverlap3}{}
#'  \item{overlap3}{}
#'  \item{nonoverlap10}{}
#'  \item{overlap10}{}
#'  }
"simdata"

#' Gene expression datasets for benchmarking classification performance
#'
#' Microarray and RNAseq datasets from biclustlib, originally published by
#' de Souto MC, Costa IG, de Araujo DS, Ludermir TB, Schliep A in 2008.
#'
#' @format A list of 13 datasets where each element contains:
#' \describe{
#'   \item{data}{a dataframe of expression values, unmodified from biclustlib}
#'   \item{labels}{a binary matrix where rows are pre-classified samples and
#'     columns are named cancer categories}
#' }
#' @source \url{https://github.com/padilha/biclustlib/tree/master/}
"cancer_benchmark"

#' Gene expression datasets for benchmarking detection of biological signal
#'
#' A dataset containing the microarray datasets in biclustlib originally
#' published by Jaskowiak PA, Campello RJ, Costa Filho IG (2013).
#'
#' @format A list of 17 data.frames where columns are samples and rows are
#'   genes, with ENSEMBL names
#'
#' @source \url{https://github.com/padilha/biclustlib/tree/master/}
#' 
"yeast_benchmark"
