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
          source("R/fileUI.R", local = TRUE)$value,
          #### Bicluster tabpanel ####
          source("R/biclusterUI.R", local = TRUE)$value,
          tabPanel(
            "Optimize",
            sidebarLayout(
              sidebarPanel(
                actionButton("bcvButton", "Perform BCV", ),
                actionButton("bcvAndBiclusterButton", "Perform BCV and bicluster", ),
                renderUI("goBicluster"),
                width = 3),
              mainPanel(
                p(id = "bcvtable", ""),
                plotOutput("bcvPlot"),
                width = 9),
              position = "left")
          ), 
          tabPanel(
            "Functional Annotation",
            sidebarLayout(
              sidebarPanel(
                actionButton("go", "Search for GO enrichment", disabled = TRUE),
                width = 3),
              mainPanel(
                tabsetPanel(
                  tabPanel("Biclusters",  plotOutput("goSigPlot")),
                  tabPanel("Terms", fluidRow(DT::DTOutput("goTermTable")))),
                # column(width = 12, 
                # p(id = "hello", "HERE@S SOME TEXT"),
                # fluidRow(DT::DTOutput("goTermTable")),
                       # style = "overflow-y: scroll"),
                width = 9), position = "left")),
          id = "navbarpage")
        )
    },
    
    #### SERVER ##################################################################
    server = source("R/server.R", local = TRUE)$value,
    
    # launch App in a browser
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})