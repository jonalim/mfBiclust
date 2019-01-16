tabPanel(
  "Data",
  tags$head(tags$style(type="text/css", # Enables auto width for rendered abundance
                       "#abundance{max-width: 100%; width: 100%; height: 100vh !important;}")
            ),
  #### sidebar ####
  sidebarLayout(
    sidebarPanel(
      h3("Upload options"),
      fileInput("input_df",label = NULL,
                accept=c('text/csv','text/comma-separated-values,text/plain',
                         '.csv')),

      checkboxInput("row_names", "Row names?", FALSE),
      checkboxInput("header", "Header?", FALSE),
      div(numericInput("skiplines", "Skip additional lines", 0, min = 0),
          style = "font-size:80%;"),
      div(textInput("sepchar", "Sep. (empty = whitespace)", ""),
          style = "font-size:80%;"),
      div(textInput("quotechar", "Quote", ""),
          style = "font-size:80%;"),
      div(textInput("decchar", "Decimal", "."),
          style = "font-size:80%;"),
      width = 3
    ),
    #### main panel and description column ####
    
    mainPanel( 
      tabsetPanel(
        tabPanel("Table",
                 tags$style(type = "text/css", "#dt {height: calc(90vh - 80px) !important;}"),
                 sidebarLayout(
                   sidebarPanel(
                     h3("Formatting"),
                     checkboxInput("transpose", "Transpose (samples must be in columns)"),
                     div(textAreaInput("customRowNames", "Custom feature names",
                                       placeholder = "Line-separated names;\nENSEMBL IDs required for GO analysis",
                                       resize = "vertical"),
                         style = "font-size:80%;"),
                     div(textAreaInput("customColNames", "Custom sample names", resize = "vertical"),
                         style = "font-size:80%;"),
                     actionButton(inputId = "postUploadUpdate", label = "Update"),
                     width = 3
                   ),
                   mainPanel(DT::DTOutput("dt"), width = 9),
                   position = "right"
                 )),
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
                width = 3
              ),
              mainPanel(
                fluidRow(uiOutput("uiabundance", width = "100%")),
                width = 9),
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
    )
  ),
  position = "left", fluid = FALSE, class = "outerTabPanel"
)
