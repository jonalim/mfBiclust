% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bicluster.R
\name{nipals_pca}
\alias{nipals_pca}
\title{Principal component dimensionality reduction using NIPALS}
\usage{
nipals_pca(A, k, cleanParam = 0, verbose = TRUE, ...)
}
\arguments{
\item{A}{the matrix to factorize}

\item{k}{the number of factors to compute}

\item{cleanParam}{passed to \code{\link{clean}()}}

\item{verbose}{report recursive calls and all values of \code{cleanParam}}

\item{...}{Additional parameters will be passed to
\code{\link[nipals]{nipals}}.}
}
\value{
a list containing
  \describe{
  \item{m}{the data matrix after any cleaning}
  \item{genericFit}{a \code{\link{genericFit-class}} object}
  \item{indexRemaining}{a list of the row and column indexes remaining after
  cleaning}
  }
}
\description{
Factorizes matrix \code{A} as the product of score and loading matrices
respectively truncated to \code{k} rows and \code{k} columns. Uses the
Nonlinear Iterative Partial Least Squares algorithm to compute principal
components in the presence of missing matrix elements.
}
\details{
If NIPALS fails, this function will recursively call itself with decreasing
values of \code{cleanParam} until NIPALS succeeds.
}
