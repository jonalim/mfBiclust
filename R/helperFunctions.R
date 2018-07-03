#### capitalize ####
capitalize <- Vectorize(function(s) {
  # Use this function whenever names will be displayed in plots or gui
  if (s == "als-nmf") { "ALS-NMF" }
  if (s == "nipals-pca") { "NIPALS-PCA" }
  if (s == "svd-pca") { "SVD-PCA" }
  else { switch(s, snmf = "SNMF", pca = "PCA", otsu = "Otsu", 
                paste0(toupper(substring(s, 1,1)), substring(s, 2))) }
})

#' Create annotations dataframe for heatmaps
#'
#' Checks phenoLabels and biclustLabels for validity and then joins the
#' corresponding data extracted from the provided BiclusterExperiment.
#'
#' Annotation tracks are converted to named dataframe columns.
#'
#' @param x a BiclusterExperiment
#' @param names the names of entities annotated- required
#' @param strategy a strategy in x. Required whenever length(biclustLabels) > 0
#' @param phenoLabels any phenotype labels in x
#' @param biclustLabels any bicluster labels in x
createAnnots <-
  function(x,
           names,
           strategy = "",
           phenoLabels = c(),
           biclustLabels = c()) {
    annots <- NA
    # Process phenotype labels
    phenoLabels <- validatePhenoNames(phenoLabels, x)
    if (length(phenoLabels) > 0) {
      phData <-
        as.data.frame(Biobase::pData(Biobase::phenoData(x))[, phenoLabels])
      colnames(phData) <- phenoLabels
    }
    
    # Process bicluster labels
    if (length(strategy) > 0 && length(biclustLabels) == 0) {
      warning(
        paste0(
          "Since bicluster annotations were not ",
          "requested, the \"strategy\" argument will be ignored."
        )
      )
    } else if (length(biclustLabels) > 0) {
      validateStratName(strategy, x) # if no strategy, stop
      bcs <- getStrat(x, strategy) # get BiclusterStrategy
      biclustLabels <- validateBiclustNames(biclustLabels, bcs)
      predData <- as.data.frame(pred(bcs)[, biclustLabels])
      
      # bicluster annotations should be discrete
      predData[] <- as.data.frame(lapply(predData, as.factor))
      colnames(predData) <- biclustLabels
    }
    
    # Concatenate phenotype and bicluster labels if both requested
    annots <-
      if (length(phenoLabels) > 0 && length(biclustLabels) > 0) {
        cbind(phData, predData)
      } else if (length(phenoLabels) > 0) {
        phData
      }
    else if (length(biclustLabels) > 0) {
      predData
    }
    
    # pheatmap throws cryptic error without rownames
    if (inherits(annots, "data.frame")) {
      row.names(annots) <- names
    }
    annots
  }

duplicable <- function() {
  if (!exists(".Random.seed", mode="numeric")) sample(NA)
  oldSeed <- .Random.seed
  set.seed(12345)
  oldSeed
}

#' Error bars
#' 
#' Function adapted from James Holland Jones, 2009.
#' http://monkeysuncle.stanford.edu/?p=485
#' 
#' @param x the result of a call to barplot()
#' @param y the ys provided to barplot()
#' @param upper vector of lengths of upper error bars
#' @param lower vector of lengths of lower error bars; by default, same as 
#'   upper
#' @export
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  suppressWarnings(arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...))
}

filter.biclust <- function(RowxBicluster, BiclusterxCol, max = NULL, 
                           overlap = 0.25) {
  if(ncol(RowxBicluster) != nrow(BiclusterxCol)) { # validate
    stop(paste0("RowxBicluster must have the same number of columns as",
                "BiclusterxCol"))
  }
  
  k <- ncol(RowxBicluster)
  if(k == 1 || k == 0) {
    # no filtering needed
    list(RowxBicluster = RowxBicluster, 
         BiclusterxCol = BiclusterxCol)
  } else {
    
    # Create lists of rows and columns contained in biclusters
    BiclusterRows <- apply(RowxBicluster, MARGIN = 2, which)
    BiclusterCols <- apply(BiclusterxCol, MARGIN = 1, which)
    
    chosen <- rep(FALSE, each = k)
    pool <- rep(TRUE, each = k)
    sizes <- sapply(seq_len(k), function(biclus) {
      length(BiclusterRows[[biclus]]) * length(BiclusterCols[[biclus]])
    })
    names(sizes) <- seq_len(k)
    
    # exclude empty biclusters and whole-dataset biclusters
    pool[sizes == 0 | sizes == nrow(RowxBicluster) * ncol(BiclusterxCol)] <- FALSE
    
    # For each bicluster, calculate its overlap with all other biclusters
    overlaps <- lapply(seq_len(k), function(biclus1) {
      sapply(seq_len(k), function(biclus2) {
        intersection1 <- base::intersect(BiclusterRows[biclus1], BiclusterRows[biclus2])
        intersection1 <- if(length(intersection1) > 0) length(intersection1[[1]])
        else 0
        intersection2 <- base::intersect(BiclusterCols[biclus1], BiclusterCols[biclus2])
        intersection2 <- if(length(intersection2) > 0) length(intersection2[[1]])
        else 0
        intersection <- intersection1 * intersection2
        overlap <- intersection / sizes[biclus1]
      })
    })
    
    while(all(sum(chosen) < max) && sum(pool) > 0) {
      chooseMe <- as.numeric(names(which.max(sizes[which(pool)])))
      chosen[chooseMe] <- TRUE
      # Remove from the pool any biclusters heavily overlapping with chooseMe
      pool[overlaps[[chooseMe]] > overlap] <- FALSE
      pool[chooseMe] <- FALSE
    }
    list(RowxBicluster = RowxBicluster[, chosen], 
         BiclusterxCol = BiclusterxCol[chosen, ])
  }
}

is.wholenumber <-
  function(x, tol = sqrt(.Machine$double.eps)) {
    abs(x - round(x)) < tol
  }
# 
# autoNipals <- function(m, k, cleanParam = 0) {
#   cleanRes <- clean(m, cleanParam, index = TRUE)
#   mClean <- cleanRes$obj
#   indexRem <- cleanRes$indexRemaining
#   
#   tryCatch({
#     list(m = mClean, nipals_pca(mClean, k), indexRemaining = indexRem)
#   }, error = function(e) {
#     if(grepl(pattern = paste0("replacement has length zero"), x = e)) {
#       cleanParam <- cleanParam + (1 - cleanParam) / 2
#       message(paste("Too many NA in the data. Cleaning with maxNAs at", 
#                     cleanParam))
#       # pass the original m so indexRemaining is valid for the user's matrix
#       autoNipals(m, k, cleanParam)
#     } else { stop(e) }
#   })
# }

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

# spectral may find over k biclusters, but only k will be returned
# minSize can be used to force biclusters to be a certain fraction of the smaller matrix dimension
spectral <- function(m, k, minSize = NULL, reps = 1) {
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  minx <- if(is.null(minSize)) 2 else floor(min(nrow(m), ncol(m)) * minSize)
  number <- 0 # save the biclustering solution with the most clusters
  
  # JNL try to find the lowest value of withinVar that yields enough biclusters.
  # 10 * nrow(m) is an arbitrary cutoff. Note that this is most effective when
  # height is the smaller matrix dimension.
  v <- nrow(m) 
  best <- NULL
  while(number < k && v <= 10L * nrow(m)) {
    bc <- biclust::biclust(m, method = biclust::BCSpectral(), normalization = "log",
                           withinVar= v, minr = minx, minc = minx)
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
  
  scores <- sapply(seq_len(k), function(i, biclusters) {
    bicluster <- biclusters[[i]]
    s <- rep(0, nrow(m))
    s[bicluster$Rows] <- 1
    s
  }, biclusters = biclusters)
  
  loadings <- do.call(rbind, lapply(seq_len(k), function(i, biclusters) {
    bicluster <- biclusters[[i]]
    l <- rep(0, ncol(m))
    l[bicluster$Cols] <- 1
    l
  }, biclusters = biclusters))
  
  assign(".Random.seed", oldSeed, envir=globalenv()) # reset the R environment
  
  new("genericFit", fit = new("genericFactorization",
                              W = scores, H = loadings), method = "spectral")
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
  
  scores <- sapply(seq_len(k), function(i, biclusters) {
    bicluster <- biclusters[[i]]
    s <- rep(0, nrow(m))
    s[bicluster$Rows] <- 1
    s
  }, biclusters = biclust::biclusternumber(best))
  
  loadings <- do.call(rbind, lapply(seq_len(k), function(i, biclusters) {
    bicluster <- biclusters[[i]]
    l <- rep(0, ncol(m))
    l[bicluster$Cols] <- 1
    l
  }, biclusters = biclust::biclusternumber(best)))
  
  new("genericFit", fit = new("genericFactorization",
                              W = scores, H = loadings), method = "plaid")
  # write function to filter?
  # use biclust::isoverlapp res$Overlapping will be TRUe or FALSe
  # if FALSE, get the two largest biclusters
  # else, get the alrgest, then get the largest bicluster with <0.25 overlap
}

#' Wrapper for prcomp
#'
#' Returns a \code{\link{genericFit}} object.
#'
#' @param m the target matrix
#' @param k the number of principal components
#' @export
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
#' @export
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

#### threshold ####
#' Apply threshold to a score or loading matrix
#'
#' Returns a binary matrix of the same size as \code{m} where all elements over
#' the threshold are 1.
#'
#' If th is a vector, the first element of th will be used as threshold for the
#' first col/row in m, etc.
setGeneric("threshold", signature = c("m", "th"), function(m, th, ...) {standardGeneric("threshold")})

#' @export
setMethod("threshold", c(m = "matrix", th = "numeric"), function(m, MARGIN = 2, th) {
  if(length(th) == 1) {
    mat <- matrix(TRUE, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
    mat[m < th] <- FALSE
    mat
  }
  else {
    if(MARGIN == 1) {
      if(length(th) != nrow(m)) { stop("Length of th must equal nrow(m).") }
      mat <- do.call(rbind, lapply(seq_len(nrow(m)), function(row) {
        m[row, ] > th[row]
      }))
    } else {
      if(length(th) != ncol(m)) { stop("Length of th must equal ncol(m)") }
      mat <- do.call(cbind, lapply(seq_len(ncol(m)), function(col) {
        m[, col] > th[col]
      }))
    }
    colnames(mat) <- colnames(m)
    rownames(mat) <- rownames(m)
    mat
  }
}
)

setMethod("threshold", c(m = "matrix", th = "matrix"), function(m, MARGIN = 2, th) {
  threshold(m, MARGIN, as.numeric(th))
}
)

validateBiclustNames <- function(biclustNames, bcs) {
  if (length(biclustNames) > 0) {
    validNames <- names(bcs)
    if (any(!biclustNames %in% validNames)) {
      warning(
        paste0(
          "Requested bicluster labels ",
          paste(setdiff(biclustNames, validNames), sep = ", "),
          " do not match any names in the \"pred\" slot of ",
          "the named BiclusterStrategy. Call \"?names\"",
          "to see how to view bicluster names."
        )
      )
      biclustNames <- intersect(biclustNames, validNames)
    }
    biclustNames
  } else {
    biclustNames
  }
}


validateK <- function(k) {
  if (!is.atomic(k)) {
    warning(
      "Argument \"k\" contains multiple elements. Attempting to use the
      first \"k\" provided."
    )
    k <- k[1]
  }
  
  if (!is.numeric(k)) {
    warning("Attempting to coerce \"k\" to numeric. Please inspect output
            carefully!")
    k <- as.numeric(k)
    if (is.na(k)) {
      stop("Couldn't coerce \"k\" to numeric.")
    }
  }
  
  if (!all(is.wholenumber(k))) {
    warning("Rounding \"k\" to whole number. Please inspect output carefully!")
    k <- round(k)
  }
  
  if (!(k > 0)) {
    stop("Argument \"k\" must be positive.")
  }
  
  k
}

validateKM <- function(k, m = NULL, method) {
  k <- validateK(k)
  if(method == "als-nmf" || "method" == "svd-pca" || method == "nipals-pca") {
    if (k >= min(nrow(m), ncol(m))) {
      warning(paste("Initializing k to the size of the smaller matrix",
                    "dimension."))
      return(min(nrow(m), ncol(m)))
    }
  }
  if(method == "plaid" || method == "spectral") k
  k
}

validateKs <- function(ks) {
  unlist(lapply(ks, function(k)
    validateK(k)))
}

validatePhenoNames <- function(phenoNames, bce) {
  if (length(phenoNames) > 0) {
    # If provided, phenoLabels must be in phenoData labels.
    validNames <- colnames(Biobase::phenoData(bce))
    if (any(!phenoNames %in% validNames)) {
      warning(
        paste0(
          "Requested phenoLabels ",
          paste(setdiff(phenoNames, validNames), sep = ", "),
          " are not in the phenoData slot of object x."
        )
      )
      phenoNames <- intersect(phenoNames, validNames)
    }
    phenoNames
  } else {
    phenoNames
  }
}


validateStratName <- function(stratName, bce) {
  validNames <- names(bce)
  if (length(stratName) > 1) {
    stop("More than one BiclusterStrategy name provided.")
  }
  if (length(stratName) == 0) {
    stop("BiclusterStrategy name was expected but not provided.")
  }
  if (!stratName %in% validNames) {
    stop(
      paste0(
        "Argument \"strategy\" is not the name of any",
        "BiclusterStrategy object in the provided ",
        "BiclusterExperiment."
      )
    )
  }
}

###% Adapted from NMF v0.21.0 written by Renaud Gaujoux, Cathal Seoighe. (2018)
###% https://cran.r-project.org/web/packages/NMF/
###% http://renozao.github.io/NMF
#' @importFrom NMF .fcnnls
als_nmf <- function(A, x, maxIter= 100L, eta=0, beta=0.00, bi_conv=c(0, 10), eps_conv=1e-4, verbose=FALSE){
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  m = nrow(A); n = ncol(A); erravg1 = numeric();
  
  #eta=param[1]; beta=param[2]; 
  maxA=max(A); if ( eta<0 ) eta=maxA;
  eta2=eta^2;
  
  sqrteps <- sqrt(.Machine$double.eps)
  
  # bi_conv
  if( length(bi_conv) != 2 )
    stop("SNMF/", version, "::Invalid argument 'bi_conv' - value should be a 2-length numeric vector")
  wminchange=bi_conv[1]; iconv=bi_conv[2];
  
  ## VALIDITY of parameters
  # eps_conv
  if( eps_conv <= 0 )
    stop("SNMF/", version, "::Invalid argument 'eps_conv' - value should be positive")
  # wminchange
  if( wminchange < 0 )
    stop("SNMF/", version, "::Invalid argument 'bi_conv' - bi_conv[1] (i.e 'wminchange') should be non-negative")
  # iconv
  if( iconv < 0 )
    stop("SNMF/", version, "::Invalid argument 'bi_conv' - bi_conv[2] (i.e 'iconv') should be non-negative")
  # beta
  # if( beta <=0 )
  #   stop("SNMF/", version, "::Invalid argument 'beta' - value should be positive")
  # ##
  
  # initialize random W if no starting point is given
  if( is.numeric(x) ){
    # rank is given by x
    k <- x
    message('# NOTE: Initialise W internally (runif)')
    W <- matrix(runif(m*k), m,k);	
    
    x <- NULL
  } else if( is.nmf(x) ){
    # rank is the number of basis components in x
    k <- nbasis(x)
    # seed the method (depends on the version to run)
    start <- t(coef(x))
    # check compatibility of the starting point with the target matrix
    if( any(dim(start) != c(m,k)) )
      stop("SNMF/", version, " - Invalid initialization - incompatible dimensions [expected: ", paste(c(m,k), collapse=' x '),", got: ", paste(dim(start), collapse=' x '), " ]")	
    # use the supplied starting point
    W <- start
  }else{
    stop("SNMF/", version, ' - Invalid argument `x`: must be a single numeric or an NMF model [', class(x), ']')
  }
  
  if ( verbose )
    cat(sprintf("--\nAlgorithm: SNMF/%s\nParameters: k=%d eta=%.4e beta (for sparse H)=%.4e wminchange=%d iconv=%d\n",
                version, k,eta,beta,wminchange,iconv));
  
  idxWold=rep(0, m); idxHold=rep(0, n); inc=0;
  
  # check validity of seed
  if( any(NAs <- is.na(W)) )
    stop("SNMF/", version, "::Invalid initialization - NAs found in the ", if(version=='R') 'basis (W)' else 'coefficient (H)' , " matrix [", sum(NAs), " NAs / ", length(NAs), " entries]")
  
  # normalize columns of W
  W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );
  Wold <- W
  Hold <- matrix(runif(k*n), k,n);	
  dnormOld <- maxA
  
  I_k=diag(eta, k); betavec=rep(sqrt(beta), k); nrestart=0;
  i <- 0L
  while( i < maxIter){
    i <- i + 1L
    
    # min_h ||[[W; 1 ... 1]*H  - [A; 0 ... 0]||, s.t. H>=0, for given A and W.
    res = .fcnnls(rbind(W, betavec), rbind(A, rep(0, n)))
    H = res[[1]]
    
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
      W=matrix(runif(m*k), m,k);
      W= apply(W, 2, function(x) x / sqrt(sum(x^2)) );  # normalize columns of W	
      next;
    }
    
    # min_w ||[H'; I_k]*W' - [A'; 0]||, s.t. W>=0, for given A and H. 
    res = .fcnnls(rbind(t(H), I_k), rbind(t(A), matrix(0, k,m))); 
    Wt = res[[1]]
    W= t(Wt);		
    
    # track the error (not computed unless tracking option is enabled in x)
    if( !is.null(x) ) 
      x <- trackError(x, .snmf.objective(A, W, H, eta, beta), niter=i)
    
    # test convergence every 5 iterations OR if the base average error has not been computed yet
    if ( (i %% 5==0)  || (length(erravg1)==0) ){
      
      #### Convergence test adapted from:####
      ###% M.W. Berry et al. (2007), "Algorithms and Applications for Approximate
      ###% Nonnegative Matrix Factorization," Computational Statistics and Data
      ###% Analysis, vol. 52, no. 1, pp. 155-173.
      dnorm = sqrt(sum((A - W %*% H)^2) / length(A))
      dw = max(abs(W - Wold) / (sqrteps + max(abs(Wold))))
      dh = max(abs(H - Hold) / (sqrteps + max(abs(Hold))))
      delta = max(dw, dh)
      if(delta <= eps_conv) break
      else if(dnormOld - dnorm <= eps_conv * max(1, dnormOld)) break
      # end adapted
      
      # # indice of maximum for each row of W
      # idxW = max.col(W)
      # # indice of maximum for each column of H
      # idxH = max.col(t(H))
      # changedW=sum(idxW != idxWold); changedH=sum(idxH != idxHold);
      # if ( (changedW<=wminchange) && (changedH==0) ) inc=inc+1
      # else inc=0
      # 
      # resmat=pmin(H, crossprod(W) %*% H - t(W) %*% A + matrix(beta, k , k) %*% H); resvec=as.numeric(resmat);
      # resmat=pmin(W, W %*% tcrossprod(H) - A %*% t(H) + eta2 * W); resvec=c(resvec, as.numeric(resmat));
      # conv=sum(abs(resvec)); #L1-norm      
      # convnum=sum(abs(resvec)>0);
      # erravg=conv/convnum;
      # # compute base average error if necessary
      # if ( length(erravg1)==0 )
      #   erravg1=erravg;
      
      if ( verbose && (i %% 100==0) ){ # prints number of changing elements
        cat("Track:\tIter\tdeltaMaxChange\tdNorm\n")
        cat(sprintf("\t%d\t%f\t%f\n",
                    i,delta, dnorm))
      }
      
      # #print(list(inc=inc, iconv=iconv, erravg=erravg, eps_conv=eps_conv, erravg1=erravg1))
      # if ( (inc>=iconv) && (erravg<=eps_conv*erravg1) ) break;
      # idxWold=idxW; idxHold=idxH; 
      Hold <- H
      Wold <- W
    }
    
  }
  
  if( verbose ) cat("--\n")
  
  # force to compute last error if not already done
  if( !is.null(x) ) 
    x <- trackError(x, .snmf.objective(A, W, H, eta, beta), niter=i, force=TRUE)	
  
  # transpose and reswap the roles
  if( !is.null(x) ){ 
    .basis(x) <- t(H)
    .coef(x) <- t(W)
    
    # set number of iterations performed
    niter(x) <- i
    
    return(x)	
  }else{
    res <- new("genericFactorization", W= W, H = H)
    res <- new("genericFit", fit = res, method = "als-nmf")
    return(invisible(res))
  }
}
