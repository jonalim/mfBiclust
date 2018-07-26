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
      width = 9),
    position = "left")
)
