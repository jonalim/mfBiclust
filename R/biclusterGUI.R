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
  
  userBce <- obj
  
  shinyApp(
    ui = fluidPage(
      shinyjs::useShinyjs(),
      #### UI ##################################################################
      tags$head(tags$style(type="text/css", # Enables auto width for rendered image1
      "#image1 img {max-width: 100%; width: 100%; height: 100%}")),
      navbarPage(
        "mfBiclust UI",
        tabPanel( #### Data tabpanel ####
                  "Data",
                  source("R/ui.R", local = TRUE)$value),
        tabPanel( #### Bicluster tabpanel ####
                  "Bicluster",
                  sidebarLayout(
                    sidebarPanel(
                      selectInput("algo", label = "Biclustering algorithm",
                                  choices = c("ALS-NMF", "SVD-PCA", "NIPALS-PCA",
                                              "SNMF", "Plaid", "Spectral")
                      ),
                      uiOutput("kSlider"),
                      actionButton("bicluster", "Run", )
                    ),
                    mainPanel(
                      tabsetPanel(
                        tabPanel("Summary", {
                          # filled interactive heatmap
                          # In a imageOutput, passing values for click, dblclick, hover, or brush
                          # will enable those interactions.
                          imageOutput("image1",
                                      # Equivalent to: click = clickOpts(id = "image_click")
                                      click = "image_click",
                                      dblclick = dblclickOpts(
                                        id = "image_dblclick"
                                      ),
                                      hover = hoverOpts(
                                        id = "image_hover"
                                      ),
                                      brush = brushOpts(
                                        id = "image_brush"
                                      )
                          )
                        }),
                        tabPanel("Inspect samples", {
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("selectCluster"),
                              checkboxInput("scoreReorder", "Reorder"),
                              # uiOutput("annotPicker"),
                              checkboxInput("sampNames", "Sample names"),
                              width = 2),
                            mainPanel(
                              uiOutput("uiScoreHeatmap"),
                              plotOutput("score_threshold", width = "100%"),
                              width = 10),
                            position = "right")
                          # Score-thresholded heatmap (try empty heatmap with annotations?)
                        }),
                        tabPanel("Inspect features", {})
                      )
                    )
                  )
        ), collapsible = FALSE, fluid = FALSE, id = "navbar"
      )
    ),
    
    #### SERVER ##################################################################
    server = source("R/server.R", local = TRUE)$value,
    
    # launch App in a browser
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})