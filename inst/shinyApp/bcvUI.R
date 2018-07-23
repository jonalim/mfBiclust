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
)
