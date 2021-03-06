% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{strategies}
\alias{strategies}
\title{Results of biclustering runs on a single dataset}
\usage{
strategies(x)
}
\arguments{
\item{x}{A \code{BiclusterExperiment-class} object}
}
\value{
A list containing zero or more \code{\link{BiclusterStrategy-class}}
  objects
}
\description{
Returns the contents of \code{\link{BiclusterExperiment-class}@strategies}.
}
\examples{
bce <- BiclusterExperiment(yeast_benchmark[[1]])
bce <- addStrat(bce, k = 2, method = "als-nmf")
bce <- addStrat(bce, k = 3, method = "als-nmf")
strategies(bce)
}
