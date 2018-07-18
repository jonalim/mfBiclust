tabPanel(
  "Functional Annotation",
  sidebarLayout(
    sidebarPanel(
      actionButton("go", "Test for GO enrichment", disabled = TRUE),
      conditionalPanel(condition = "input.goTab == 'Terms' ||
                                 input.goTab == 'Genes'",
                       uiOutput("goBicluster")),
      conditionalPanel(condition = "input.goTab == 'Genes'",
                       selectInput("goTerm", label = "Select GO ID:", choices = list(),
                                   multiple = FALSE)),
      width = 3),
    mainPanel(
      tabsetPanel(id = "goTab",
                  tabPanel("Biclusters",  plotOutput("goSigPlot")),
                  tabPanel("Terms", 
                           fluidRow(DT::DTOutput("goTermTable"),
                                    style = "height:500px; overflow-y: scroll"),
                           actionButton("goTabGenes",
                                        label = "Inspect genes")),
                  tabPanel("Genes", 
                           column(6,
                                  h6("Matched genes in bicluster:"),
                                  verbatimTextOutput("goBiclusterGenes")),
                           column(6,
                                  h6("Matched genes in dataset:"),
                                  verbatimTextOutput("goUniverseGenes")))),
      width = 9), 
    position = "left"))
