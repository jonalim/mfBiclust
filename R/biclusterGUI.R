#' @include BiclusterExperiment.R
NULL

# Please remember to use clusterProfiler for testing GO enrichment

#' Explore a BiclusterExperiment
#'
#' Opens a GUI to visualize analysis results contained in a BiclusterExperiment
#' object
#'
#' @export
setGeneric("biclusterGUI", signature = "obj", function(obj = NULL) {
  standardGeneric("biclusterGUI")
})

# clusters must be a named list of matrices
#' @describeIn biclusterGUI Open GUI for a BiclusterExperiment
#' @import shiny
setMethod("biclusterGUI", definition = function(obj) {
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
          tabPanel("Data", 
                   tags$head(tags$style(type="text/css", # Enables auto width for rendered image1
                                        "#image1 img {max-width: 100%; width: 100%; height: 100%}")),
                   source("R/fileUI.R", local = TRUE)$value),
          #### Bicluster tabpanel ####
          source("R/biclusterUI.R", local = TRUE)$value,
          tabPanel(
            "Optimize",
            sidebarLayout(
              sidebarPanel(
                actionButton("bcvButton", "Perform BCV", ),
                actionButton("bcvAndBiclusterButton", "Perform BCV and bicluster", ),
                width = 3),
              mainPanel(
                p(id = "bcvtable", ""),
                plotOutput("bcvPlot"),
                width = 9),
              position = "left")
          ), id = "navbarpage"
        )
      )
    },
    
    #### SERVER ##################################################################
    server = source("R/server.R", local = TRUE)$value,
    
    # launch App in a browser
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})