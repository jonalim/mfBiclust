#' Biclustering algorithms
#'
#' The \code{method} argument of \code{\link{addStrat}()} provides access to
#' functions from packages \code{\link{[NMF]}} and \code{\link{[biclust]}}.
#' Initially \code{addStrat()} uses the default parameters provided by the
#' respective developers, then progressively relaxes certain parameters if
#' needed to return the desired number of biclusters. Method-specific parameters
#' besides those automatically manipulated can be provided by name to
#' \code{addStrat()}. Please see respective documentation in other packages for
#' method-specific parameters.
#' 
#' \describe{ 
#' \item{als_nmf}{The alternating-least-squares non-negative matrix
#' factorization algorithm. The original algorithm by Paatero and Tapper (1994)
#' is run \code{reps} times, then the result that most accurately factorized the
#' matrix \code{A} is returned. \code{als_nmf} is fast for a small number of
#' bicluster, but running time rapidly increases with \code{k}.}
#' \item{svd_pca}{The singular vector decomposition algorithm. Each principal
#' component is interpreted as the degree of membership in a single bicluster.
#' The resulting score matrix is thresholded to binarize bicluster membership.
#' \code{svd_pca} is the fastest provided algorithm.}
#' \item{nipals_pca}{An iterative PCA algorithm that may tolerate missing data.
#'   Slower than \code{svd_pca}, but still faster than the other algorithms.}
#' \item{plaid}{The plaid algorithm}
#' \item{snmf}{Sparse non-negative matrix factorization. The core algorithm is
#' Paatero and Tapper's alternating least-squares (1994), and the
#' sparsity-inducing regularization factors were introduced by Kim and Park
#' (2007).}
#' \item{spectral}{The Spectral algorithm. For smaller matrices, the parameter
#'   \code{numberOfEigenvalues} is automatically set to enable finding the
#'   number of biclusters requested. The parameter \code{withinVar} is allowed
#'   to increase up to 10 times the smaller matrix dimension.}
#' }
#' 
#' @param A a numeric matrix to bicluster
#' @param k the number of biclusters to report
#' @param duplicable causes biclutsering to be non-stochastic
#' @param verbose report details such as convergence and the underlying function
#'   calls
#'   
#' @seealso \code{\link[NMF]{nmfAlgorithm.SNMF_R}}
#' @seealso \code{\link[biclust]{BCSpectral}}
#' @seealso \code{\link[biclust]{BCPlaid}}
#' @seealso \code{\link[nipals]{nipals}}
#' @name bicluster
NULL

#' @importFrom NMF .fcnnls
als_nmf <- function(A, k, reps = 4L, maxIter= 100L, eta=0L, beta=0.00, 
                    eps_conv = 1e-7, duplicable = TRUE, verbose=TRUE, ...){
  ###% Adapted from NMF v0.21.0 written by Renaud Gaujoux, Cathal Seoighe.
  ###% (2018)
  ###% https://cran.r-project.org/web/packages/NMF/
  ###% http://renozao.github.io/NMF
  if(duplicable) {
    oldSeed <- duplicable("biclus") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  } else { reps <- 1 }
  
  m = nrow(A); n = ncol(A); erravg1 = numeric();
  
  # to do with sparsity constraints
  maxA=max(A); if ( eta<0 ) eta=maxA;
  eta2=eta^2;
  
  ## VALIDITY of parameters
  # eps_conv
  if( eps_conv <= 0 )
    stop("SNMF/", version, "::Invalid argument 'eps_conv' - value should be positive")
  # beta
  if( beta < 0 )
    stop("SNMF/", version, "::Invalid argument 'beta' - value should be positive")
  
  solutions <- lapply(seq_len(reps), function(i) {
    W <- NULL # these are the results of each replicate
    H <- NULL
    residNorm <- max(A)
    # initialize random W if no starting point is given
    
    # rank is given by k
    W <- matrix(runif(m*k), m,k)
    
    idxWold=rep(0, m); idxHold=rep(0, n); inc=0;
    
    # normalize columns of W
    frob <- apply(W, 2, function(x) sqrt(sum(x ^ 2)) )
    W <- sweep(W, 2, frob, FUN = "/")
    
    Wold <- W
    Hold <- matrix(runif(k*n), k,n);	
    residNormOld <- maxA
    I_k=diag(eta, k); betavec=rep(sqrt(beta), k); nrestart=0;
    restart <- function() {
      # re-initialize random W
      idxWold <<- rep(0, m); idxHold <<- rep(0, n); inc <<- 0; 
      erravg1 <<- numeric();# re-initialize base average error
      W <- matrix(runif(m*k), m,k);
      
      # normalize columns of W	
      frob <- apply(W, 2, function(x) sqrt(sum(x ^ 2)) )
      W <<- sweep(W, 2, frob, FUN = "/")
    }
    
    i <- 0L
    while( i < maxIter){
      i <- i + 1L
      
      # min_h ||[[W; 1 ... 1]*H  - [A; 0 ... 0]||, s.t. H>=0, for given A and W.
      res <- try(.fcnnls(rbind(W, betavec), rbind(A, rep(0, n))), silent = TRUE)
      if(inherits(res, "try-error")) {
        nrestart <- nrestart+1;
        if ( nrestart >= 10 ){
          warning("NMF::snmf - Factorization failed too many times [Computation stopped after the 9th restart]");
          break;
        }
        restart()
        next
      }
      H <- res[[1]]
      
      if ( any(rowSums(H)==0) ){
        nrestart <- nrestart+1;
        if ( nrestart >= 10 ){
          warning("NMF::snmf - Factorization failed too many times [Computation stopped after the 9th restart]");
          break;
        }
        restart()
        next
      }
      
      # min_w ||[H'; I_k]*W' - [A'; 0]||, s.t. W>=0, for given A and H.
      res = try(.fcnnls(rbind(t(H), I_k), rbind(t(A), matrix(0, k,m))),
                silent = TRUE)
      if(inherits(res, "try-error")) {
        nrestart <- nrestart+1;
        if ( nrestart >= 10 ){
          warning("NMF::snmf - Factorization failed too many times [Computation stopped after the 9th restart]");
          break;
        }
        restart()
        next
      }
      
      Wt = res[[1]]
      W <- t(Wt);		
      
      #### Convergence test adapted from:####
      ###% M.W. Berry et al. (2007), "Algorithms and Applications for Approximate
      ###% Nonnegative Matrix Factorization," Computational Statistics and Data
      ###% Analysis, vol. 52, no. 1, pp. 155-173.
      residNorm <- sum((A - W %*% H) ^ 2) / length(A)
      if(abs(residNormOld - residNorm) <= eps_conv && verbose) {
        if(verbose) message("Converged!")
        break
      }
      residNormOld <- residNorm
      # end adapted
      
      Hold <- H
      Wold <- W
      
      # every 5 iterations
      if ( (i %% 5==0)  || (length(erravg1)==0) ){
        if ( verbose && (i %% 5==0) ){ # prints number of changing elements
          cat("Track:\tIter\tNorm\tdelta\n")
          cat(sprintf("\t%d\t%f\t%f\n",
                      i,residNorm, residNormOld - residNorm))
        }
      }
    }
    return(list(obj = residNorm, W = W, H = H))
  })
  
  solution.best <- which.min(unlist(sapply(solutions, function(x) x$obj)))
  
  W <- solutions[[solution.best]]$W
  H <- solutions[[solution.best]]$H
  
  res <- new("genericFactorization", W= W, H = H)
  res <- new("genericFit", fit = res, method = "als-nmf")
  return(res)
}

nipals_pca_nocatch <- function(A, k, duplicable = TRUE, center = FALSE, ...) {
  if(duplicable) {
    oldSeed <- duplicable("biclus") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  }
  
  np <- tryCatch({
    nipals::nipals(x = A, ncomp = k, center = center, scale = FALSE, 
                   tol = 1e-6)
  },
  error = function(e) {
    if(grepl(pattern = paste0("replacement has length zero"), x = e)) {
      stop(paste("NIPALS was not able to run. Try running clean() on your",
                 "input matrix or BiclusterExperiment object."))
    } else {
      stop(e)
    }
  })
  new("genericFit", fit = new("genericFactorization",
                              W = np$scores,
                              H = t(np$loadings)),
      method = "nipals-pca")
}

nipals_pca <- function(A, k, cleanParam = 0, center = FALSE,
                                 duplicable = FALSE, verbose = TRUE, ...) {
  if(duplicable) {
    oldSeed <- duplicable("biclus") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  }
  cleanRes <- clean(A, cleanParam, dimsRemain = TRUE)
  mClean <- cleanRes$obj
  indexRem <- cleanRes$dimsRemain
  
  tryCatch({
    list(m = mClean, genericFit = nipals_pca_nocatch(mClean, k, center = center),
         indexRemaining = indexRem)
  }, error = function(e) {
    if(grepl(pattern = paste0("replacement has length zero"), x = e)) {
      cleanParam <- cleanParam + log10(2 - i)
      if(verbose) {
        message(paste("Too many NA in the data. Cleaning with maxNAs at",
                    cleanParam))
      }
      # pass the original m so indexRemaining is valid for the user's matrix
      nipals_pca_helper(A, k, cleanParam, center, duplicable)
    } else { stop(e) }
  })
}

plaid <- function(A, k, duplicable = TRUE, verbose = TRUE, ...) {
  if(duplicable) {
    oldSeed <- duplicable("biclus") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  }
  
  args <- list(...)
  shuffleMax <- args$shuffle
  args <- c(cluster = args$cluster, fit.model = args$fit.model, 
            background = args$background, 
            background.layer = args$background.layer,
            background.df = args$background.df,
            backfit = args$backfit, iter.startup = args$iter.startup, 
            iter.layer = args$iter.layer, verbose = args$verbose)
  if(is.null(shuffleMax)) shuffleMax <- 10
  number <- 0
  release <- 0.7
  shuff <- 3
  best <- NULL
  while(shuff <= shuffleMax && number < k) {
    while(number < k && release > 0) {
      dummy <- capture.output({
        bc <- do.call(biclust::biclust, 
                      c(list(x = A, method = biclust::BCPlaid(),
                             row.release = release, col.release = release,
                             max.layers = k, shuffle = shuff), args))
      })
      if(bc@Number > number) {
        number <- bc@Number
        best <- bc
      }
      release <- max(0.1, release - 0.1)
      # first release decrements from 0.7 to 0.1
    }
    # then if necessary, shuffle increments from 3 to 10
    shuff <- shuff + 1
  }
  
  if(k > number) {
    k <- number
    warning(paste("Plaid could only find", k, "biclusters"))
  }
  if(verbose) {
    cat(paste("method =", class(bc@Parameters$Call$method), "\n"))
    cat(paste("row.release =", bc@Parameters$Call$row.release, "\n"))
    cat(paste("col.release =", bc@Parameters$Call$col.release, "\n"))
    cat(paste("max.layers =", bc@Parameters$Call$max.layers, "\n"))
    cat(paste("shuffle =", bc@Parameters$Call$shuffle, "\n"))
  }
  scoreLoading <- if(k > 0) { 
    biclusters <- biclust::biclusternumber(best)
    biclusterNumber2scoreLoading(biclusters, A, k) 
  } else { list(matrix(rep(NA, nrow(A)), ncol = 1), 
                matrix(rep(NA, ncol(A)), nrow = 1)) }
  
  new("genericFit", fit = new("genericFactorization",
                              W = scoreLoading[[1]], H = scoreLoading[[2]]), method = "plaid")
}

snmf <- function(A, k, beta = 0.01, verbose = TRUE, duplicable = TRUE, ...) {
  if(duplicable) {
    oldSeed <- duplicable("biclus") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  }
  
  args <- list(...)
  args <- c(maxIter = args$maxIter, eta = args$eta,bi_conv = args$bi_conv, 
            eps_conv = args$eps_conv, .options = args$.options)
  
  res <- tryCatch(
    {
      # /r causes columns (sample membership) to be sparse. Presumably,
      # samples are more likely than features to correspond 1:1 with a bicluster
      suppressMessages(
        do.call(NMF::nmf, c(list(x = A, rank = k, method = "snmf/r", 
                                 beta = beta, rng = .Random.seed), args)))
    },
    warning = function(w) {
      if (any(suppressWarnings(
        grepl(
          "too big 'beta' value",
          w$message,
          ignore.case = TRUE,
          fixed = TRUE
        )
      ))) {
        beta <<- beta ^ 2
        res <<- snmf(A, k, beta, verbose, duplicable, ...)
      } else {
        warning(w)
      }
    },
    error = function(e) {
      stop(e)
    }
  )
  if(verbose) {
    cat(paste("method:", res@method, "\n"))
    cat(paste("parameters:\n", res@parameters, "\n"))
  }
  res@method <- "snmf"
  return(res)
}

spectral <- function(A, k, minSize = NULL, reps = 1, duplicable = TRUE, 
                     verbose = TRUE, ...) {
  # spectral may find over k biclusters, but only k will be returned
  # minSize can be used to force biclusters to be a certain fraction of the smaller matrix dimension
  if(duplicable) {
    oldSeed <- duplicable("biclus") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  }
  
  minx <- if(is.null(minSize)) 2 else floor(min(nrow(A), ncol(A)) * minSize)
  if(nrow(A) < 6 || ncol(A) < 6) {
    stop("For Spectral the minimum size of m is 6x6")
  }
  number <- 0 # save the biclustering solution with the most clusters
  
  withinVar <- list(...)$withinVar
  
  # JNL try to find the lowest value of withinVar that yields enough biclusters.
  # 10 * nrow(m) is an arbitrary cutoff designed for when rows
  # are samples. The user can, however, override withinVar. The GUI uses this
  # parameter in debug mode to limit running time.
  if(is.null(withinVar)) {
    withinVar <- min(nrow(A), ncol(A))
  }
  
  # the number of eigenvalues to consider. This can limit number of biclusters
  # found, so increase from the default if necessary.
  e <- 3
  if(!nrow(A) / minx * ncol(A) / minx * e >= k) {
    e <- ceiling(k / (nrow(A) / minx * ncol(A) / minx))
  }
  
  best <- NULL
  while(number < k && withinVar <= 10L * min(nrow(A), ncol(A))) {
    bc <- do.call(biclust::biclust,
                  list(x = A, method = biclust::BCSpectral(),
                       normalization = "log", withinVar= withinVar, minr = minx,
                       minc = minx, numberOfEigenvalues = e))
    if(bc@Number > number) {
      number <- bc@Number
      best <- bc
    }
    withinVar <- withinVar + nrow(A)
  }
  if(k > number) {
    k <- number
    warning(paste("Spectral could only find", k, "biclusters"))
  }
  
  if(verbose) {
    cat(paste("method =", class(bc@Parameters$Call$method), "\n"))
    cat(paste("normalization =", bc@Parameters$Call$normalization, "\n"))
    cat(paste("withinVar =", bc@Parameters$Call$withinVar, "\n"))
    cat(paste("minr =", bc@Parameters$Call$minr, "\n"))
    cat(paste("minc =", bc@Parameters$Call$minc, "\n"))
    cat(paste("numberOfEigenvalues =", bc@Parameters$Call$numberOfEigenvalues, "\n"))
  }
  
  scoreLoading <- if(k > 0) { 
    biclusters <- biclust::biclusternumber(best)
    biclusterNumber2scoreLoading(biclusters, A, k) 
  } else { list(matrix(rep(NA, nrow(A)), ncol = 1), 
                matrix(rep(NA, ncol(A)), nrow = 1)) }
  
  new("genericFit", fit = new("genericFactorization", W = scoreLoading[[1]], 
                              H = scoreLoading[[2]]),
      method = "spectral")
}

svd_pca <- function(A, k, duplicable = NULL, ...) {
  prcmp <- prcomp(t(A), rank. = k, retx = TRUE, center = FALSE)
  new(
    "genericFit",
    fit = new(
      "genericFactorization",
      H = t(prcmp$x),
      W = prcmp$rotation
    ),
    method = "svd-pca"
  )
}
