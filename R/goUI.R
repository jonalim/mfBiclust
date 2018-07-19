tabPanel(
  "Functional Annotation",
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("gos", h6("Ontologies"),
                         choices = list("Biological Process" = "BP",
                                        "Molecular Function" = "MF",
                                        "Cellular Component" = "CC"),
                         selected = "BP"),
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
                           tags$style(type = "text/css", "#goTermTable {height: calc(90vh - 80px) !important;}"),
                           DT::DTOutput("goTermTable"),
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
