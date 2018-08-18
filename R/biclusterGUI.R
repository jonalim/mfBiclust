#' @include BiclusterExperiment.R
NULL

#' Explore a BiclusterExperiment
#'
#' Opens a shiny GUI to visualize analysis results contained in a
#' BiclusterExperiment object
#' 
#' @param obj a BiclusterExperiment, an ExpressionSet, or an object
#'   coercible to a matrix
#' @param dbg Runs in a debug mode that uses algorithm parameters increasing
#'   sacrificing accuracy for speed.
#'
#' @example R/examples/addStrat-biclusterGUI.R
#' @export
#' @name biclusterGUI
#' @import shiny
setGeneric("biclusterGUI", signature = "obj",
           function(obj = NULL, ...) 
             standardGeneric("biclusterGUI"))

setMethod("biclusterGUI", c(obj = "BiclusterExperiment"), 
          definition = function(obj, ...) {
  ## define UI parameters
  plotHeight <- 600
  plotHeightSmall <- 300
  
  dbg <- list(...)$dbg
  if(is.null(dbg)) { dbg <- FALSE }
  
  # what server.R is expecting, if no good BiclusterExperiment is available
  if(all(is.na(as.matrix(obj)))) { obj <- NULL }
  userBce <- obj # server.R has access to this variable
  
  # Display infinity as "Inf" in json (without this, Inf values completely
  # missing in datatables)
  tojson_args.old <- getOption("DT.TOJSON_ARGS")
  options("DT.TOJSON_ARGS" = list(na = 'string'))
  on.exit(expr = options("DT.TOJSON_ARGS" = tojson_args.old), add = TRUE)
  
  shinyApp(
    ui = source(system.file("shinyApp", "ui.R", package = "mfBiclust"),
                local = TRUE)$value,
    server = source(system.file("shinyApp", "server.R", package = "mfBiclust"),
                    local = TRUE)$value,
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})

setMethod("biclusterGUI", c(obj = "ExpressionSet"), 
          definition = function(obj, ...) {
            biclusterGUI(as(obj, "BiclusterExperiment"), ...)
          })
setMethod("biclusterGUI", c(obj = "missing"), function(...) {
  biclusterGUI(obj = BiclusterExperiment(m = matrix()), ...)
})
setMethod("biclusterGUI", c(obj = "ANY"), 
          definition = function(obj, ...) {
            biclusterGUI(BiclusterExperiment(obj), ...)
          })
