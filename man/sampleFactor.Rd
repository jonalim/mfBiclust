% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{sampleFactor}
\alias{sampleFactor}
\title{Sample matrix factor}
\usage{
sampleFactor(bcs, allBc = FALSE)
}
\arguments{
\item{bcs}{A BiclusterStrategy}

\item{allBc}{TRUE to return data from all stored biclusters}
}
\value{
a numeric matrix.
}
\description{
For a data matrix \eqn{A_{m x n}} factorized as \eqn{W_{m x k}H_{k x n}},
returns \eqn{H}.
}
\examples{
bce <- BiclusterExperiment(yeast_benchmark[[1]])
bce <- addStrat(bce, k = 2, method = "als-nmf")
bcs <- getStrat(bce, 1)
sampleFactor(bcs)

}
