sidebarLayout( #### sidebar ####
               sidebarPanel(
                 fileInput("input_df", "Import data", accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
                 checkboxInput("sampleCols", "Samples as columns"),
                 checkboxInput("row_names", "Row names?", FALSE),
                 checkboxInput("header", "Header?", FALSE),
                 numericInput("skiplines", "Skip additional lines", 0),
                 textInput("sepchar", "Sep. (empty = whitespace)", ""),
                 textInput("quotechar", "Quote", ""),
                 textInput("decchar", "Decimal", "."),
                 width = 3
               ),
               mainPanel( #### main panel and description column ####
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
                                         selectInput("distType", "Formula",
                                                     choices = c("Euclidean distance",
                                                                 "1 - Pearson")), 
                                         checkboxInput("distReorder", "Reorder"),
                                         width = 3),
                                       mainPanel(plotOutput("distance", width = "100%"),
                                                 width = 9),
                                       position = "right")
                            )
                          ),
                          width = 9
               ), fluid = FALSE
)