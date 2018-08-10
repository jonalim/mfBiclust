#### Bicluster Tab Panel for the Shiny UI

tabPanel(
  "Bicluster",
  tags$head(
    tags$style(HTML(
      "label { font-size: 120% ; }
      .navbar { font-size: 120% ; }
      .navbar-brand { font-size: 120% ; }"
    ))
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("algo", label = "Biclustering algorithm",
                  choices = c("ALS-NMF" = "als-nmf", "SVD-PCA" = "svd-pca", 
                              "NIPALS-PCA" = "nipals-pca", "SNMF" = "snmf",
                              "Plaid" = "plaid", "Spectral" = "spectral"),
                  multiple = FALSE
      ),
      sliderInput("k", label = "Number of biclusters:", min = 1, max = 2,
                  value = 1, step = 1),
      actionButton("bicluster", "Run"),
      uiOutput("downloadBceButton"),
      uiOutput("downloadBceAllButton"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Summary",
          tags$style(
            type = "text/css",
            paste0("#bcHighlights {image-rendering: pixelated !important;
                   image-rendering: optimizeSpeed;
          image-rendering: -moz-crisp-edges;
          image-rendering: -o-crisp-edges;
          image-rendering: -webkit-optimize-contrast;
          image-rendering: optimize-contrast;")),
          # Try removing optimize contrast or optimizeSpeed if the heatmap
          # doesn't look right
          
          # filled interactive heatmap
          # In a imageOutput, passing values for click, dblclick, hover, or brush
          # will enable those interactions.
          imageOutput("bcHighlights",
                      # Equivalent to: click = clickOpts(id = "bcHighlights_click")
                      click = "bcHighlights_click",
                      dblclick = dblclickOpts(
                        id = "bcHighlights_dblclick"
                      ),
                      hover = hoverOpts(
                        id = "bcHighlights_hover"
                      ),
                      brush = brushOpts(
                        id = "bcHighlights_brush"
                      )
          )
        ),
        tabPanel("Samples",
                 column(9,
                        plotOutput("sampleHeatmap", width = "100%"), 
                        plotOutput("samplePlot", width = "100%")),
                 column(3,
                        uiOutput("sampleBicluster"),
                        checkboxInput("sampleReorder", "Reorder"),
                        checkboxInput("biclusterFeatNames", "Sample names"))
        ),
        tabPanel("Features",
                 column(9,
                        plotOutput("featureHeatmap", width = "100%"), 
                        plotOutput("featurePlot", width = "100%")),
                 column(3,
                        uiOutput("featureBicluster"),
                        checkboxInput("featureReorder", "Reorder"),
                        # uiOutput("annotPicker"),
                        checkboxInput("biclusterSampNames", "Feature names"),
                 uiOutput("biclusterGeneListLabel"),
                 fluidRow(verbatimTextOutput("biclusterGeneList"), 
                          style = "height:500px; overflow-y: scroll"))
        )
      ),
      width = 9, 
      
      position = "left")
  ), class = "outerTabPanel"
)
