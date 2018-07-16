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
                choices = c("ALS-NMF", "SVD-PCA", "NIPALS-PCA",
                            "SNMF", "Plaid", "Spectral")
    ),
    uiOutput("kSlider"),
    actionButton("bicluster", "Run", ),
    width = 3
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
      tabPanel("Samples",
               column(9,
                      plotOutput("scoreHeatmap", width = "100%"), 
                      plotOutput("scorePlot", width = "100%")
               ), column(
                 3,
                 uiOutput("scoreBicluster"),
                 checkboxInput("scoreReorder", "Reorder"),
                 # uiOutput("annotPicker"),
                 checkboxInput("sampNames", "Sample names"))
               # Score-thresholded heatmap (try empty heatmap with annotations?)
      ),
      tabPanel("Features", {
        fluidRow(
          column(9,
                 plotOutput("loadingHeatmap", width = "100%"), 
                 plotOutput("plot_biomarkers", width = "100%")),
          column(3,
                 uiOutput("loadingBicluster"),
                 checkboxInput("loadingReorder", "Reorder"),
                 checkboxInput("featNames", "Feature names"),
                 uiOutput("biclusterGeneListLabel"),
                 fluidRow(verbatimTextOutput("biclusterGeneList"), 
                          style = "height:500px; overflow-y: scroll"))
        )
        # uiOutput("loadingPanel")
      })),
      width = 9
    ),
    position = "left")
  )