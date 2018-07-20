# RETURN list of named vectors of BH-adjusted p-values
# Each list element should be named after a bicluster
testFE <- function(bce, strategy, orgDb = "org.Sc.sgd.db", go = c("BP", "MF"), duplicable = TRUE) {
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
  
  # Change this to use method dispatch
  bcs <- if(inherits(strategy, "BiclusterStrategy")) strategy else getStrat(bce, strategy)
  
  if(duplicable) { duplicable("testFE") }
  
  # thresholded loading matrix
  biclusterxCol <- threshold(loading(bcs), MARGIN = 1, th = loadingThresh(bcs))
  
  geneLists <- lapply(seq_len(nrow(biclusterxCol)), function(biclust) {
    rownames(bce)[biclusterxCol[biclust, ] == 1]
  })
  names(geneLists) <- names(strategy)
  
  universe <- rownames(bce) # these must be ENSEMBL
  
  # Returns a list of results for each bicluster: a list with one element per
  # ontology
  if(length(go) > 0) {
    mFun <- function(geneList, name, universe, ontology, fun, orgDb) {
      # test each requested GO
      tryCatch(
        {BiocParallel::bplapply(go, fun, 
              geneList = geneList, name = name, universe = universe,
              orgDb = orgDb)
        },
        error = function(e) { # fallback if BiocParallel failes
          lapply(go, fun, 
                 geneList = geneList, name = name, universe = universe,
                 orgDb = orgDb)
        }
      )
    }
    return(# Test every bicluster for enrichment
      tryCatch(
        {BiocParallel::bpmapply(mFun,
              geneLists, names(geneLists),
              MoreArgs = list(universe = universe, ontology = go,
                              fun = hyperGGO, orgDb = orgDb),
              SIMPLIFY = FALSE, USE.NAMES = TRUE)
        },
        error = function(e) {
          mapply(mFun,
                 geneLists, names(geneLists),
                 MoreArgs = list(universe = universe, ontology = go,
                                 fun = hyperGGO, orgDb = orgDb),
                 SIMPLIFY = FALSE, USE.NAMES = TRUE)
        }
      )
    )
  }
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