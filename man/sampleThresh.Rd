% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{sampleThresh}
\alias{sampleThresh}
\title{Sample thresholds}
\usage{
sampleThresh(bcs, allBc = FALSE)
}
\arguments{
\item{bcs}{A BiclusterStrategy}

\item{allBc}{TRUE to return data from all stored biclusters}
}
\value{
a numeric vector of length \eqn{k}
}
\description{
Returns the thresholds that were applied to each row of the factor matrix
\eqn{H_{k x n}} to determine the samples in each bicluster.
}
\examples{
bce <- BiclusterExperiment(yeast_benchmark[[1]])
bce <- addStrat(bce, k = 2, method = "als-nmf")
bcs <- getStrat(bce, 1)
sampleThresh(bcs)
}
