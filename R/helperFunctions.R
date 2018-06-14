

#' Capitalize names and abbreviations
#'
#' Use this function whenever names will be displayed in plots or gui
#' 
capitalize <- Vectorize(function(s) {
  if (s == "als-nmf") { "ALS-NMF" }
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

is.wholenumber <-
  function(x, tol = sqrt(.Machine$double.eps)) {
    abs(x - round(x)) < tol
  }

#' Wrapper for prcomp
#'
#' Returns a \code{\link{genericFit}} object.
#'
#' @param m the target matrix
#' @param k the number of principal components
#' @export
pcaWrapper <- function(m, k) {
  prcmp <- prcomp(m, rank. = k, retx = TRUE)
  new(
    "genericFit",
    fit = new(
      "genericFactorization",
      W = prcmp$x,
      H = t(prcmp$rotation)
    ),
    method = "pca"
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

# Apply threshold to a score or loading matrix
setGeneric("threshold", signature = "m", function(m, ...) {standardGeneric("threshold")})
#' @export
setMethod("threshold", c(m = "matrix"), function(m, th) {
  mat <- matrix(TRUE, nrow = nrow(m), ncol = ncol(m))
  mat[m < th] <- FALSE
  mat
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

validateKM <- function(k, m = NULL) {
  k <- validateK(k)
  if (k >= nrow(m) || k >= ncol(m)) {
    stop(
      paste0(
        "Number of biclusters must be smaller than both dimensions of",
        "the assay data matrix."
      )
    )
  }
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

# Adapted from NMF v0.21.0 written by Renaud Gaujoux, Cathal Seoighe. (2018)
# https://cran.r-project.org/web/packages/NMF/
# http://renozao.github.io/NMF
#'
#' @importFrom NMF .fcnnls
als_nmf <- function(A, x, maxIter= 100L, eta=0, beta=0.00, bi_conv=c(0, 10), eps_conv=1e-4, verbose=FALSE){
  #nmfsh_comb <- function(A, k, param, verbose=FALSE, bi_conv=c(0, 10), eps_conv=1e-4, version=c('R', 'L')){
  
  # # depending on the version: 
  # # in version L: A is transposed while W and H are swapped and transposed
  # version <- match.arg(version)
  # if( version == 'L' ) A <- t(A) 
  #if( missing(param) ) param <- c(-1, 0.01)
  
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
      
       #### Algorithm adapted from:####
       #     M.W. Berry et al. (2007), "Algorithms and Applications for Approximate
       #     Nonnegative Matrix Factorization," Computational Statistics and Data
       #     Analysis, vol. 52, no. 1, pp. 155-173.
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
