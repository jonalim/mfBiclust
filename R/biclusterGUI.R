#' @include BiclusterExperiment.R
NULL

#' Explore a BiclusterExperiment
#'
#' Opens a shiny GUI to analyze and visualize bicluster analyses.
#'
#' @param obj a BiclusterExperiment, an ExpressionSet, or an object
#'   coercible to a matrix. If missing, the user must import data from a comma-
#'   or space-delimited file using the GUI.
#' @param ... Other parameters. \code{dbg = TRUE} Runs the GUI in a debug mode
#'   that uses algorithm parameters tweaked to sacrifice accuracy for speed.
#'
#' @return A Shiny app object
#'
#' @examples
#' \dontrun{biclusterGUI()
#' biclusterGUI(yeast_benchmark[[1]])
#' biclusterGUI(BiclusterExperiment(yeast_benchmark[[1]]))
#' }
#'
#' @export
#' @name biclusterGUI
#' @import shiny shinythemes
#' @importFrom shinyjs useShinyjs click enable disable show hide html runjs
setGeneric("biclusterGUI", signature = "obj",
           function(obj = NULL, ...)
             standardGeneric("biclusterGUI"))
#' @describeIn biclusterGUI Default method
#' @importFrom DT DTOutput renderDT datatable
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
#' @describeIn biclusterGUI Coerces \code{obj} to a
#'   \code{\link{BiclusterExperiment-class}} object and runs GUI
setMethod("biclusterGUI", c(obj = "ExpressionSet"),
          definition = function(obj, ...) {
            biclusterGUI(as(obj, "BiclusterExperiment"), ...)
          })
#' @describeIn biclusterGUI Runs GUI without pre-loading a dataset. Use the
#'   "Data" tab to load data.
setMethod("biclusterGUI", c(obj = "missing"), function(...) {
  biclusterGUI(obj = BiclusterExperiment(m = matrix()), ...)
})
#' @describeIn biclusterGUI Attempts to encapsulate \code{obj} in a
#'   \code{\link{BiclusterExperiment-class}} object and pass it to the GUI
setMethod("biclusterGUI", c(obj = "ANY"),
          definition = function(obj, ...) {
            biclusterGUI(BiclusterExperiment(obj), ...)
          })
