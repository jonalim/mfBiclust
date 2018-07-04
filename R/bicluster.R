###% Adapted from NMF v0.21.0 written by Renaud Gaujoux, Cathal Seoighe. (2018)
###% https://cran.r-project.org/web/packages/NMF/
###% http://renozao.github.io/NMF
#' @importFrom NMF .fcnnls
als_nmf <- function(A, x, rep = 4, maxIter= 100L, eta=0, beta=0.00, 
                    eps_conv = sqrt(.Machine$double.eps), verbose=FALSE){
  # oldSeed <- duplicable() # do not modify the R global environment
  # on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  m = nrow(A); n = ncol(A); erravg1 = numeric();
  
  # to do with sparsity constraints
  maxA=max(A); if ( eta<0 ) eta=maxA;
  eta2=eta^2;
  
  sqrteps <- sqrt(.Machine$double.eps)
  
  ## VALIDITY of parameters
  # eps_conv
  if( eps_conv <= 0 )
    stop("SNMF/", version, "::Invalid argument 'eps_conv' - value should be positive")
  # beta
  if( beta <=0 )
    stop("SNMF/", version, "::Invalid argument 'beta' - value should be positive")
  
  solutions <- lapply(seq_len(4), function(x) {
    W <- NULL # these are the results of each replicate
    H <- NULL
    residNorm <- max(A)
    # initialize random W if no starting point is given
    
      # rank is given by x
      k <- x
      message('# NOTE: Initialise W internally (runif)')
      W <- matrix(runif(m*k), m,k)
      
      x <- NULL
    
    idxWold=rep(0, m); idxHold=rep(0, n); inc=0;
    
    # check validity of seed
    if( any(NAs <- is.na(W)) )
      stop("SNMF/", version, "::Invalid initialization - NAs found in the ", if(version=='R') 'basis (W)' else 'coefficient (H)' , " matrix [", sum(NAs), " NAs / ", length(NAs), " entries]")
    
    # normalize columns of W
    frob <- apply(W, 2, function(x) sqrt(sum(x ^ 2)) )
    W <- sweep(W, 2, frob, FUN = "/")
    
    Wold <- W
    Hold <- matrix(runif(k*n), k,n);	
    residNormOld <- maxA
    I_k=diag(eta, k); betavec=rep(sqrt(beta), k); nrestart=0;
    i <- 0L
    while( i < maxIter){
      i <- i + 1L
      
      # min_h ||[[W; 1 ... 1]*H  - [A; 0 ... 0]||, s.t. H>=0, for given A and W.
      res = .fcnnls(rbind(W, betavec), rbind(A, rep(0, n)))
      H <- res[[1]]

      if ( any(rowSums(H)==0) ){
        if( verbose ) cat(sprintf("iter%d: 0 row in H eta=%.4e restart!\n",i,eta));
        nrestart=nrestart+1;
        if ( nrestart >= 10 ){
          warning("NMF::snmf - Too many restarts due to too big 'beta' value [Computation stopped after the 9th restart]");
          break;
        }
        
        # re-initialize random W
        idxWold=rep(0, m); idxHold=rep(0, n); inc=0; 
        erravg1 <- numeric();# re-initialize base average error
        W <- matrix(runif(m*k), m,k);
        
        # JNL normalize columns of W	
        frob <- apply(W, 2, function(x) sqrt(sum(x ^ 2)) )
        W <- sweep(W, 2, frob, FUN = "/")
        
        next
      }
      
      # min_w ||[H'; I_k]*W' - [A'; 0]||, s.t. W>=0, for given A and H. 
      res = .fcnnls(rbind(t(H), I_k), rbind(t(A), matrix(0, k,m))); 
      Wt = res[[1]]
      W <- t(Wt);		
      
      # track the error (not computed unless tracking option is enabled in x)
      if( !is.null(x) ) 
        x <- trackError(x, .snmf.objective(A, W, H, eta, beta), niter=i)
      
      #### Convergence test adapted from:####
      ###% M.W. Berry et al. (2007), "Algorithms and Applications for Approximate
      ###% Nonnegative Matrix Factorization," Computational Statistics and Data
      ###% Analysis, vol. 52, no. 1, pp. 155-173.
      residNorm <- sum((A - W %*% H) ^ 2)
      if(abs(residNormOld - residNorm) <= eps_conv) {
        message("Converged!")
        break
      }
      residNormOld <- residNorm
      # end adapted
      
      Hold <- H
      Wold <- W
      
      # every 5 iterations
      if ( (i %% 5==0)  || (length(erravg1)==0) ){
        if ( verbose && (i %% 5==0) ){ # prints number of changing elements
          cat("Track:\tIter\tdeltaMaxChange\tdNorm\n")
          cat(sprintf("\t%d\t%f\t%f\n",
                      i,0, residNorm))
        }
      }
    }
    return(list(obj = residNorm, W = W, H = H))
  })

  solution.best <- which.min(unlist(sapply(solutions, function(x) {
    
    x$obj
    })))
  W <- solutions[[solution.best]]$W
  H <- solutions[[solution.best]]$H
  
  res <- new("genericFactorization", W= W, H = H)
  res <- new("genericFit", fit = res, method = "als-nmf")
  return(invisible(res))
}
nipals_pca <- function(m, k, reps = 1) {
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  np <- tryCatch({
    nipals::nipals(x = m, ncomp = k, center = FALSE, scale = FALSE, 
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

plaid <- function(m, k) {
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  number <- 0
  release <- 0.7
  best <- NULL
  while(number < k && release > 0) {
    dummy <- capture.output({
      bc <- biclust::biclust(m, method = biclust::BCPlaid(), 
                             row.release = release, col.release = release,
                             max.layers = k)
    })
    if(bc@Number > number) {
      number <- bc@Number
      best <- bc
    }
    release <- release - 0.1
  }
  
  if(k > number) {
    k <- number
    warning(paste("Plaid could only find", k, "biclusters"))
  }
  biclusters <- biclust::biclusternumber(best)
  
  scoreLoading <- if(k > 0) { biclusterNumber2scoreLoading(biclusters, m, k) }
  else { list(matrix(rep(NA, nrow(m)), ncol = 1), 
              matrix(rep(NA, ncol(m)), nrow = 1)) }
  
  new("genericFit", fit = new("genericFactorization",
                              W = scoreLoading[[1]], H = scoreLoading[[2]]), method = "plaid")
}

#' Wrapper for NMF::nmf(method = "snmf/l")
#'
#' This wrapper explicitly enables use of the non-regularized alternating least
#' squares algorithm devised by Paatero and Tapper (1994), with optional
#' parameters providing access to coefficients of the regularization factors
#' introduced by Kim and Park (2007)
#'
#' An \code{NMFfit} object
#'
#' @param m the target matrix
#' @param k the size of the reduced dimension
#' @param beta the starting beta
snmf <- function(m, k, beta = 0.01, verbose = FALSE) {
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  tryCatch(
    suppressMessages(res <-
                       NMF::nmf(
                         m, k, method = "snmf/l", beta = beta,
                         verbose = verbose
                       )),
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
        message(paste0("Decreased beta (sparsity parameter) to ", beta))
        res <<- snmfWrapper(m, k, beta)
      } else {
        warning(w)
      }
    },
    error = function(e) {
      stop(e)
    },
    finally = function() {
      res
    }
  )
}

# spectral may find over k biclusters, but only k will be returned
# minSize can be used to force biclusters to be a certain fraction of the smaller matrix dimension
spectral <- function(m, k, minSize = NULL, reps = 1) {
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  minx <- if(is.null(minSize)) 2 else floor(min(nrow(m), ncol(m)) * minSize)
  if(nrow(m) < 6 || ncol(m) < 6) {
    stop("For Spectral the minimum size of m is 6x6")
  }
  number <- 0 # save the biclustering solution with the most clusters
  
  # JNL try to find the lowest value of withinVar that yields enough biclusters.
  # 10 * nrow(m) is an arbitrary cutoff. Note that this is most effective when
  # height is the smaller matrix dimension.
  v <- nrow(m)
  
  # the number of eigenvalues to consider. This can limit number of biclusters
  # found, so increase from the default if necessary.
  e <- 3
  if(!nrow(m) / minx * ncol(m) / minx * e >= k) {
    e <- ceiling(k / (nrow(m) / minx * ncol(m) / minx))
  }
  best <- NULL
  while(number < k && v <= 10L * nrow(m)) {
    bc <- biclust::biclust(m, method = biclust::BCSpectral(), normalization = "log",
                           withinVar= v, minr = minx, minc = minx,
                           numberOfEigenvalues = e)
    if(bc@Number > number) {
      number <- bc@Number
      best <- bc
    }
    v <- v + nrow(m)
  }
  if(k > number) {
    k <- number
    warning(paste("Spectral could only find", k, "biclusters"))
  }
  
  biclusters <- biclust::biclusternumber(best)
  scoreLoading <- if(k > 0) { biclusterNumber2scoreLoading(biclusters, m, k) }
  else { list(matrix(rep(NA, nrow(m)), ncol = 1), 
              matrix(rep(NA, ncol(m)), nrow = 1)) }
  
  new("genericFit", fit = new("genericFactorization",
                              W = scoreLoading[[1]], H = scoreLoading[[2]]), method = "spectral")
}

#' Wrapper for prcomp
#'
#' Returns a \code{\link{genericFit}} object.
#'
#' @param m the target matrix
#' @param k the number of principal components
svd_pca <- function(m, k) {
  prcmp <- prcomp(m, rank. = k, retx = TRUE, center = FALSE)
  new(
    "genericFit",
    fit = new(
      "genericFactorization",
      W = prcmp$x,
      H = t(prcmp$rotation)
    ),
    method = "svd-pca"
  )
}
