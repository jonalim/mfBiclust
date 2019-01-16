#' Simulated bicluster datasets
#'
#' Simulated matrices containing biclusters of various models. Generated using
#' \code{\link{genSimData}()}.
#'
#' @format A list of 4 matrices named:
#' \describe{
#'  \item{fiveSquareConstant}{No background, sixteen non-overlapping 5x5
#'  biclusters. Elements in biclusters have value 5; all other elements are 0.}
#'  \item{fiveSquareRowShift}{No background, two 5x5 biclusters causing
#'  row-shift effects}
#'  \item{fiveSquareOverlapped}{Three 5x5 biclusters, where the first and second
#'  biclusters overlap and thesecond and third biclusters overlap. Both overlap
#'  regions are 4 rows by 2 columns.}
#'  \item{plaid}{A 10x10 plaid bicluster}
#'  }
"simdata"

#' Gene expression datasets for benchmarking classification performance
#'
#' Pre-classified microarray and RNAseq datasets from human cancer samples,
#' originally published by de Souto MC, Costa IG, de Araujo DS, Ludermir TB,
#' Schliep A in 2008.
#'
#' Redistributed from biclustlib: A Python library of biclustering algorithms
#' and evaluation measures. Copyright (C) 2017 Victor Alexandre Padilha.
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
#' Yeast microarray datasets originally published by Jaskowiak PA, Campello RJ,
#' Costa Filho IG (2013). All rownames are ENSEMBL gene names.
#'
#' Redistributed from biclustlib: A Python library of biclustering algorithms
#' and evaluation measures. Copyright (C) 2017 Victor Alexandre Padilha.
#'
#' @format A list of 17 matrices where columns are samples and rows are
#'   genes, with ENSEMBL names
#'
#' @source \url{https://github.com/padilha/biclustlib/tree/master/}
#'
"yeast_benchmark"
