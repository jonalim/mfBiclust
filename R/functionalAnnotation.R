# RETURN list of named vectors of BH-adjusted p-values
# Each list element should be named after a bicluster
testFE <- function(bce, strategy, orgDb = "org.Sc.sgd.db", go = c("BP", "MF", "CC"), 
                   parallel = FALSE, duplicable = TRUE) {
  if(length(strategy) != 1) {
    stop(paste("Provide exactly one BiclusterStrategy object, name or index. Run",
               "names(strategies(bce)) to see BiclusterStrategy objects."))
  }
  
  # Look for the gene-to-GO database
  if(!suppressPackageStartupMessages(requireNamespace(orgDb, quietly = TRUE))) {
    stop(paste0("Package ", orgDb, " needed. Please install it. If it is ",
                "present on BioConductor, the recommended method is ",
               "BiocInstaller::biocLite(\"", orgDb, "\")"))
  }
  
  bcs <- if(inherits(strategy, "BiclusterStrategy")) { strategy
    } else getStrat(bce, strategy)
  
  go <- match.arg(go, several.ok = TRUE)

  lappl <- lapply
  mappl <- mapply
  if(parallel) {
    if(!requireNamespace("BiocParallel", quietly = TRUE)) {
      warning("BiocParallel must be on the library path to use parallelization")
    } else {
      lappl <- BiocParallel::bplapply
      mappl <- BiocParallel::bpmapply
    }
  }
  
  if(duplicable) { duplicable("testFE") }
  
  # thresholded loading matrix
  biclusterxCol <- clusteredFeatures(bcs)
  
  geneLists <- lapply(seq_len(ncol(clusteredFeatures(bcs))), function(biclust) {
    rownames(bce)[biclusterxCol[, biclust] == 1]
  })
  names(geneLists) <- bcNames(strategy)
  universe <- rownames(bce) # these must be ENSEMBL
  # Returns a list of results for each bicluster: a list with one element per
  # ontology
  mFun <- function(geneList, name, universe, ontology, fun, orgDb) {
    # test each requested GO
    lappl(go, fun, 
          geneList = geneList, name = name, universe = universe,
          orgDb = orgDb)
  }
  return(# Test every bicluster for enrichment
    mappl(mFun,
          geneLists, names(geneLists),
          MoreArgs = list(universe = universe, ontology = go,
                          fun = hyperGGO, orgDb = orgDb),
          SIMPLIFY = FALSE, USE.NAMES = TRUE)
  )
}

hyperGGO <- function(ontology = c("MF", "BP", "CC"), geneList, name, universe,
                     orgDb) {
  message(paste("Testing", name, "for Gene Ontology term enrichment")
  )
  ontology <- match.arg(ontology)
  return(suppressPackageStartupMessages(
    GOstats::hyperGTest(new("GOHyperGParams",
                            geneIds = geneList,
                            universeGeneIds = universe,
                            annotation = orgDb,
                            ontology = ontology,
                            pvalueCutoff = 1,
                            testDirection = "over", 
                            conditional = TRUE))
  ))
}