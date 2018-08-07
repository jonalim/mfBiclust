tabPanel(
  "Optimize",
  sidebarLayout(
    sidebarPanel(
      uiOutput("bcvButton"),
      uiOutput("bcvAndBiclusterButton"),
      br(),
      p("Bi-cross-validation determines the number of biclusters that",
        "minimizes the residuals when factorizing your dataset via SVD-PCA.",
        "In practice, bi-cross-validation is also roughly applicable to", 
        "ALS-NMF."),
      width = 3),
    mainPanel(
      p(id = "bcvtable", ""),
      plotOutput("bcvPlot"),
      width = 9),
    position = "left")
)
