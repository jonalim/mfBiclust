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
          imageOutput("bcHighlights", height = "600px", width = "100%",
                      # Equivalent to: click = clickOpts(id = "bcHighlights_click")
                      click = "bcHighlights_click",
                      dblclick = dblclickOpts(
                        id = "bcHighlights_dblclick"
                      ),
                      hover = hoverOpts(
                        id = "bcHighlights_hover"
                      ),
                      brush = brushOpts(
                        id = "bcHighlights_brush", delay = 10000, clip = TRUE
                      )
          )
        ),
        tabPanel("Samples",
                 column(9,
                        plotOutput("sampleHeatmap", width = "100%"),
                        HTML("Row <em>i</em> of the heatmap shows each sample's degree of association with the <em>i</em>th bicluster, on a unitless scale relative to other samples.
                                     Values in different bicluster-rows should not be compared.<br>
                                    <br><br>To determine bicluster membership, the Otsu algorithm is used to threshold each row.
                             The plot below shows which samples were thresholded into the selected biluster. A list of samples in the selected bicluster is printed to the right."),
                        plotOutput("samplePlot", width = "100%")),
                 column(3,
                        uiOutput("sampleBicluster"),
                        checkboxInput("sampleReorder", "Reorder"),
                        checkboxInput("biclusterFeatNames", "Sample names"),
                        uiOutput("biclusterSampleListLabel"),
                        fluidRow(verbatimTextOutput("biclusterSampleList"),
                                 style = "height:500px; overflow-y: scroll"))
        ),
        tabPanel("Features",
                 column(9,
                        plotOutput("featureHeatmap", width = "100%"),
                        HTML("Row <em>i</em> of the heatmap shows the phenotype that defines the <em>i</em>th bicluster.
                                   Values are unitless, and should not be compared between biclusters.
                                   <br><br>To determine bicluster membership, the Otsu algorithm is used to threshold each row.
                             The plot below shows which features were thresholded into the selected biluster. A list of features in the selected bicluster is printed to the right."),
                        plotOutput("featurePlot", width = "100%")),
                 column(3,
                        uiOutput("featureBicluster"),
                        checkboxInput("featureReorder", "Reorder"),
                        # uiOutput("annotPicker"),
                        checkboxInput("biclusterSampNames", "Feature names"),
                        uiOutput("biclusterFeatureListLabel"),
                        fluidRow(verbatimTextOutput("biclusterFeatureList"),
                          style = "height:500px; overflow-y: scroll"))
        )
      ),
      width = 9,

      position = "left")
  ), class = "outerTabPanel"
)
