#' @include BiclusterExperiment.R
NULL

setGeneric("biclusterGUI", signature = "obj", function(obj = NULL, ...) {
  standardGeneric("biclusterGUI")
})

#' Explore a BiclusterExperiment
#'
#' Opens a shiny GUI to visualize analysis results contained in a
#' BiclusterExperiment object
#' 
#' @param dbg Runs in a debug mode that uses algorithm parameters increasing
#'   sacrificing accuracy for speed.
#'
#' @example R/examples/addStrat-biclusterGUI.R
#' @export
#' @name biclusterGUI
#' @import shiny
setMethod("biclusterGUI", c(obj = "BiclusterExperiment"), 
          definition = function(obj, dbg = FALSE) {
  ## define UI parameters
  plotHeight <- 600
  plotHeightSmall <- 300
  userBce <- obj # server.R has access to this variable
  
  shinyApp(
    ui = source(system.file("shinyApp", "ui.R", package = "mfBiclust"),
                local = TRUE)$value,
    server = source(system.file("shinyApp", "server.R", package = "mfBiclust"),
                    local = TRUE)$value,
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})