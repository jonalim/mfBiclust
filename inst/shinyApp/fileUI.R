tabPanel(
  "Data", 
  tags$head(tags$style(type="text/css", # Enables auto width for rendered image1
                       "#image1 img {max-width: 100%; width: 100%; height: 100%}")),
  #### sidebar ####
  sidebarLayout( 
    sidebarPanel(
      fileInput("input_df", "Import data", accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
      checkboxInput("sampleCols", "Samples as columns"),
      checkboxInput("row_names", "Row names?", FALSE),
      checkboxInput("header", "Header?", FALSE),
      numericInput("skiplines", "Skip additional lines", 0, min = 0),
      textInput("sepchar", "Sep. (empty = whitespace)", ""),
      textInput("quotechar", "Quote", ""),
      textInput("decchar", "Decimal", "."),
      textAreaInput("customRowNames", "Custom row names",
                    placeholder = "Line-separated ENSEMBL IDs", resize = "vertical"),
      width = 3
    ),
    #### main panel and description column ####
    mainPanel( 
      tabsetPanel(
        tabPanel("Table", fluidRow(DT::DTOutput("dt"))),
        tabPanel(
          "Heatmap",
          fluidRow(
            sidebarLayout(
              sidebarPanel(
                h4("Format"),
                selectInput("logBase", "Log-transform",
                            choices = list("None" = 0, 
                                           "log[2]" = 2, 
                                           "ln" = "e", 
                                           "log[10]" = 10)),
                checkboxInput(
                  "heatmapReorder",
                  "Reorder by distance"),
                uiOutput("annotPicker"),
                checkboxInput("sampNames", "Sample names"),
                checkboxInput("featNames", "Feature names"),
                width = 2
              ),
              mainPanel(
                fluidRow(uiOutput("uiabundance", width = "100%")),
                width = 10),
              position = "right"
            )
          )),
        tabPanel("PCA", fluidRow(plotOutput("pca"))),
        tabPanel("Sample distance", 
                 sidebarLayout(
                   sidebarPanel(
                     selectInput("sampDistType", "Formula",
                                 choices = c("Euclidean distance",
                                             "1 - Pearson")), 
                     checkboxInput("sampDistReorder", "Reorder"),
                     checkboxInput("sampDistRownames", "Show row names"),
                     width = 3),
                   mainPanel(plotOutput("sampleDistance", width = "100%"),
                             width = 9),
                   position = "right")
        )
      ),
      width = 9
    ), fluid = FALSE
  ))