#' @include BiclusterExperiment.R
NULL

# Please remember to use clusterProfiler for testing GO enrichment

#' Explore a BiclusterExperiment
#'
#' Opens a GUI to visualize analysis results contained in a BiclusterExperiment
#' object
#'
#' @export
setGeneric("biclusterGUI", signature = "obj", function(obj = NULL, ...) {
  standardGeneric("biclusterGUI")
})

# clusters must be a named list of matrices
#' @describeIn biclusterGUI Open GUI for a BiclusterExperiment
#' @import shiny
setMethod("biclusterGUI", definition = function(obj, debug = FALSE) {
  ## define UI parameters
  plotHeight <- 600
  plotHeightSmall <- 300
  library("shinythemes")
  on.exit({try(detach("package:shinythemes"))},add = TRUE)
  userBce <- obj
  
  shinyApp(
    ui = {
      tagList(
        shinyjs::useShinyjs(),
        #### UI ##################################################################
        navbarPage(
          theme = shinytheme("yeti"), inverse = TRUE, "mfBiclust UI",
          #### Data tabpanel ####
          source("R/fileUI.R", local = TRUE)$value,
          #### Bicluster tabpanel ####
          source("R/biclusterUI.R", local = TRUE)$value,
          source("R/bcvUI.R", local = TRUE)$value, 
          source("R/goUI.R", local = TRUE)$value,
          id = "navbarpage", fluid = TRUE)
        )
    },
    
    #### SERVER ##################################################################
    server = source("R/server.R", local = TRUE)$value,
    
    # launch App in a browser
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})