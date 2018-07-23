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
setMethod("biclusterGUI", definition = function(obj, dbg = FALSE) {
  ## define UI parameters
  plotHeight <- 600
  plotHeightSmall <- 300
  userBce <- obj
  
  shinyApp(
    ui = source(system.file("shinyApp", "ui.R", package = "mfBiclust"),
           local = TRUE)$value,
    server = source(system.file("shinyApp", "server.R", package = "mfBiclust"),
                    local = TRUE)$value,
    options = list(launch.browser = TRUE)
  )
  
  # shinyApp(
  #   ui = {
  #     
  #   },
  #   
  #   #### SERVER ##################################################################
  #   server = source("R/server.R", local = TRUE)$value,
  #   
  #   # launch App in a browser
  #   options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  # )
})