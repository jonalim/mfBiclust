tabPanel(
  "Optimize",
  sidebarLayout(
    sidebarPanel(
      uiOutput("bcvButton"),
      uiOutput("bcvAndBiclusterButton"),
      width = 3),
    mainPanel(
      p(id = "bcvtable", ""),
      plotOutput("bcvPlot"),
      HTML("Bi-cross-validation (BCV) is used to determine the optimal number of
            biclusters for matrix-factorization-based biclustering. Since BCV is
            a stochastic procedure, BCV is repeated until the distribution of 
            results converges. The results will be tallied while running
           BCV, and then the final results will be displayed as a bar graph.
           <br><br>Technically, the singular value decomposition (used in the 
          SVD-PCA biclustering algorithm) is being
           validated, but the result may be useful for ALS-NMF as well."),
      width = 9),
    position = "left")
)
