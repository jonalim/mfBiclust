% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{simdata}
\alias{simdata}
\title{Simulated bicluster datasets}
\format{A list of 4 matrices named:
\describe{
 \item{fiveSquareConstant}{No background, sixteen non-overlapping 5x5
 biclusters. Elements in biclusters have value 5; all other elements are 0.}
 \item{fiveSquareRowShift}{No background, two 5x5 biclusters causing
 row-shift effects}
 \item{fiveSquareOverlapped}{Three 5x5 biclusters, where the first and second
 biclusters overlap and thesecond and third biclusters overlap. Both overlap
 regions are 4 rows by 2 columns.}
 \item{plaid}{A 10x10 plaid bicluster}
 }}
\usage{
simdata
}
\description{
Simulated matrices containing biclusters of various models. Generated using
\code{\link{genSimData}()}.
}
\keyword{datasets}
