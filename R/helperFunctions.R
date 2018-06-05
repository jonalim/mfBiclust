is.wholenumber <-
  function(x, tol = .Machine$double.eps ^ 0.5) {
    abs(x - round(x)) < tol
  }

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
createAnnots <- function(x, names, strategy = "", phenoLabels = c(), biclustLabels = c()) {
  annots <- NA
  # Process phenotype labels
  phenoLabels <- validatePhenoNames(phenoLabels, x)
  if(length(phenoLabels) > 0) {
    phData <- as.data.frame(Biobase::pData(Biobase::phenoData(x))[,phenoLabels])
    colnames(phData) <- phenoLabels
  }
  
  # Process bicluster labels
  if (length(strategy) > 0 && length(biclustLabels) == 0) {
    warning(paste0("Since bicluster annotations were not ",
                   "requested, the \"strategy\" argument will be ignored.")
    )
  } else if(length(biclustLabels) > 0) {
    validateStratName(strategy, x) # if no strategy, stop
    bcs <- getStrat(x, strategy) # get BiclusterStrategy 
    biclustLabels <- validateBiclustNames(biclustLabels, bcs)
    predData <- as.data.frame(pred(bcs)[, biclustLabels])
    
    # bicluster annotations should be discrete
    predData[] <- as.data.frame(lapply(predData, as.factor))
    colnames(predData) <- biclustLabels
  }
  
  # Concatenate phenotype and bicluster labels if both requested
  annots <- if(length(phenoLabels) > 0 && length(biclustLabels) > 0) {
    cbind(phData, predData)
  } else if(length(phenoLabels) > 0) { phData } 
  else if(length(biclustLabels) > 0) { predData }
  
  # pheatmap throws cryptic error without rownames
  if(inherits(annots, "data.frame")) {
    row.names(annots) <- names
  }
  annots
}

validatePhenoNames <- function(phenoNames, bce) {
  if(length(phenoNames) > 0) {
    # If provided, phenoLabels must be in phenoData labels.
    validNames <- colnames(Biobase::phenoData(bce))
    if(any(!phenoNames %in% validNames)) {
      warning(paste0("Requested phenoLabels ",
                     paste(setdiff(phenoNames, validNames), sep = ", "),
                     " are not in the phenoData slot of object x."))
      phenoNames <- intersect(phenoNames, validNames)
    }
    phenoNames
  } else {phenoNames}
}

validateBiclustNames <- function(biclustNames, bcs) {
  if(length(biclustNames) > 0) {
    validNames <- names(bcs)
    if(any(!biclustNames %in% validNames)) {
      
      warning(paste0("Requested bicluster labels ",
                     paste(setdiff(biclustNames, validNames), sep = ", "),
                     " do not match any names in the \"pred\" slot of ",
                     "the named BiclusterStrategy. Call \"?names\"",
                     "to see how to view bicluster names."))
      biclustNames <- intersect(biclustNames, validNames)
    }
    biclustNames
  } else {biclustNames}
}

validateStratName <- function(stratName, bce) {
  validNames <- names(bce)
  if(length(stratName) > 1) {
    stop("More than one BiclusterStrategy name provided.")
  }
  if(length(stratName) == 0) {
    stop("BiclusterStrategy name was expected but not provided.")
  }
  if(!stratName %in% validNames) {stop(paste0("Argument \"strategy\" is not the name of any",
                                              "BiclusterStrategy object in the provided ",
                                              "BiclusterExperiment."))
  }
}

validateK <- function(k) {
  if (!is.atomic(k)) {
    warning("Argument \"k\" contains multiple elements. Attempting to use the
first \"k\" provided.")
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
  if(k >= nrow(m) || k >= ncol(m)) {
    stop(paste0("Number of biclusters must be smaller than both dimensions of",
                "the assay data matrix."))
  }
  k
}

validateKs <- function(ks) {
  unlist(lapply(ks, function(k) validateK(k)))
}
