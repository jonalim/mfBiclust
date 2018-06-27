#' @include BiclusterExperiment.R
NULL

#' Explore a BiclusterExperiment
#'
#' Opens a GUI to visualize analysis results contained in a BiclusterExperiment
#' object
#'
#' @export
setGeneric("biclusterGUI", signature = "obj", function(obj = NULL) {
  standardGeneric("biclusterGUI")
})

# clusters must be a named list of matrices
#' @describeIn biclusterGUI Open GUI for a BiclusterExperiment
#' @import shiny
setMethod("biclusterGUI", definition = function(obj) {
  ## define UI parameters
  plotHeight <- 600
  plotHeightSmall <- 300
  
  ## define server global variables
  user <- reactiveValues(userBce = obj)
  
  shinyApp(
    ui = fluidPage(
      shinyjs::useShinyjs(),
      #### UI ##################################################################
      tags$head(tags$style(
      )),
      navbarPage("mfBiclust UI",
        tabPanel(
          "Data",
          sidebarLayout(
            sidebarPanel(
              fileInput("input_df", "Import data", accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              checkboxInput("row_names", "Row names?", FALSE),
              checkboxInput("header", "Header?", FALSE),
              numericInput("skiplines", "Skip additional lines", 0),
              textInput("sepchar", "Sep. (empty = whitespace)", ""),
              textInput("quotechar", "Quote", ""),
              textInput("decchar", "Decimal", "."),
              width = 2
            ),
            # Data panel and description column
            mainPanel(
              tabsetPanel(
                tabPanel("Table", 
                         fluidRow(DT::DTOutput("dt"))
                ),
                tabPanel(
                  "Heatmap",
                  fluidRow(
                  sidebarLayout(
                    sidebarPanel(
                      fluidRow(
                        h4("Format"),
                        selectInput("logBase", "Log-transform",
                                    choices = list("None" = 0, 
                                                   "log[2]" = 2, 
                                                   "ln" = "e", 
                                                   "log[10]" = 10)),
                        selectInput(
                          "sampOrder",
                          "Reorder by distance",
                          choices = c(
                            "Same as input" = "input",
                            "Euclidean dist by abundance" = "distance",
                            "Euclidean dist by bicluster membership" = "cluster"
                          )
                        ),
                        uiOutput("annotPicker"),
                        checkboxInput("sampNames", "Sample names"),
                        checkboxInput("featNames", "Feature names")
                      ),
                      width = 2
                    ),
                    mainPanel(
                      fluidRow(uiOutput("uiabundance", width = "100%")),
                      width = 10),
                    position = "right"
                  )
                )),
                tabPanel("PCA",
                         fluidRow(plotOutput("pca"))
                )
              ),
              # column(
              #   9,
              #   uiOutput("inside_tabs")
              # ),
              # column(
              #   3,
              #   h4("Panel description"),
              #   htmlOutput("explanation")
              # ),
              width = 9
            )
          )
        ),
        tabPanel(
          "Bicluster",
          sidebarLayout(
            sidebarPanel(
              selectInput(inputId = "strategy", label = "Biclustering strategy"
                          # FIXME
                            ),
              mainPanel()
              )
        )
      )
    ),
    
    #### SERVER ##################################################################
    server = function(input, output, session) {
      
      dep <- reactiveValues(
        pheatmap = requireNamespace("pheatmap", quietly = TRUE),
        plaid = requireNamespace("biclust", quietly = TRUE),
        spectral = requireNamespace("biclust", quietly = TRUE)
        )
      #### DATA ####
      reactiveSeparator <- observeEvent(input$input_df, {
        # If the file input changes, try automatically changing separator to ","
        if(substr(input$input_df$name, nchar(input$input_df$name) - 2,
                  nchar(input$input_df$name)) == "csv") {
          updateTextInput(session, "sepchar", value = ",")
        }
      }
      )
      
      reactiveDF <- reactive({
        # If the user has not chosen a file yet, look for a BiclusterExperiment in the execution environment of biclusterGUI
        if(is.null(input$input_df)) {
          shinyjs::disable("row_names")
          shinyjs::disable("header")
          shinyjs::disable("skiplines")
          shinyjs::disable("sepchar")
          shinyjs::disable("quotechar")
          shinyjs::disable("decchar")
          if(exists("user") && !is.null(user$userBce)) {
            as.data.frame(t(as.matrix(user$userBce)))
          } else { validate(need(FALSE, "You may import your dataset."))
          }
        }
        # If the user has chosen a file, read it and update
        else {
          shinyjs::enable("row_names")
          shinyjs::enable("header")
          shinyjs::enable("skiplines")
          shinyjs::enable("sepchar")
          shinyjs::enable("decchar")
          shinyjs::enable("quotechar")

          # force character limits
          shinyjs::runjs("$('#sepchar').attr('maxlength', 1)")
          shinyjs::runjs("$('#quotechar').attr('maxlength', 1)")
          shinyjs::runjs("$('#decchar').attr('maxlength', 1)")

          
          rows <- if(input$row_names) 1 else NULL
          
          read.table(input$input_df$datapath, header = input$header,
                        row.names = rows, fill = TRUE, comment.char = "",
                        sep = input$sepchar, quote = input$quotechar, dec = input$decchar, skip = input$skiplines)
      }
      })

      reactiveBce <- reactive({
        # Always keep the BiclusterExperiment synchronized with (itself) or the user's file
        if(is.null(input$input_df)) {
          if(is.null(user$userBce)) {
            validate(need(FALSE, "You may import your dataset."))
          } else user$userBce
        } else {
          BiclusterExperiment(t(as.matrix(reactiveDF())))
        }
        })

      output$dt <- DT::renderDT({
        return(reactiveDF())
      }, server = FALSE, selection = "none")
      
      # helper functions
      reactiveHeatmapHeight500 <- reactive({ 
        max(500, length(input$annots) * 33 +
              sum(unlist(lapply(input$annots, function(annot) {
                length(unique(cbind(Biobase::pData(Biobase::phenoData(reactiveBce())),
                                    pred(getStrat(reactiveBce(), input$strategy))
                )[, annot]))
              }))) * 22 - 106)})
      
      observeEvent(input$pheatmap, {
        withProgress(message = "Installing pheatmap...", value = 0, {
          install.packages("pheatmap")
        })
        if(requireNamespace("pheatmap", quietly = TRUE)) {
          dep$pheatmap <- TRUE
        }
        })

      # plot abundance heatmap (original data)
      output$uiabundance <- renderUI({
        if(!dep$pheatmap) {
            actionButton("pheatmap", "Install dependency \"pheatmap\"")
          } else {
        height = reactiveHeatmapHeight500()
        plotOutput("abundance", height = height)
      }
      })
      output$abundance <- renderPlot({
        gt <- reactive_abundance()
        print(gt) # printing ensures the returned gTable is drawn to Shiny
      }, height = function() {reactiveHeatmapHeight500()})
      reactive_abundance <- reactive({
        if(!is.null(reactiveBce())) {
          
          withProgress(message = "Plotting...", value = 0, {
            set.seed(1234567)
            if(input$logBase == "e") {
              logBase <- exp(1)
            } else {
              logBase <- as.numeric(input$logBase)
            }
            phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(reactiveBce())))
            # biclustLabels <- if(length(names(reactiveBce())) > 0) {
            #   intersect(input$annots, names(getStrat(reactiveBce(), input$strategy)))
            plot(reactiveBce(), logBase = logBase, phenoLabels = phenoLabels, ordering = input$sampOrder, strategy = input$strategy, 
                 rowNames = input$sampNames, colNames = input$featNames)
          })
        }
      })
      
      # Plot of samples along first two PCs
      output$pca <- renderPlot({ 
        pca(reactiveBce())
        }
        )
      
      # Euclidean Distance between samples (applied to raw data...might be good to scale first?)
      output$distance <- renderPlot({
        gt <- reactiveDistance()
        print(gt)
      }, height = function() {reactiveHeatmapHeight500()})
      reactiveDistance <- reactive({
        validate(need(
          all(input$annot.rows %in% colnames(phenoData(reactiveBce()))),
          paste0(
            "\nPlease choose annotations from: ",
            paste(colnames(phenoData(reactiveBce())), collapse = ", ")
          )
        ))
        withProgress(message = "Plotting...", value = 0, {
          set.seed(1234567)
          phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(reactiveBce())))
          biclustLabels <- intersect(input$annots, names(getStrat(reactiveBce(), input$strategy)))
          distType <- if(input$distType == "Euclidean distance") {
            "euclidean"} else {"pearson"}
          plotDist(reactiveBce(), distType = distType, phenoLabels = phenoLabels, 
                   biclustLabels = biclustLabels, 
                   ordering = input$sampOrder, strategy = input$strategy, 
                   rowColNames = input$sampNames
          )
        })
      })

      ####


      reactiveHeatmapHeight300 <- reactive({ 
        max(300, length(input$annots) * 33 +
              sum(unlist(lapply(input$annots, function(annot) {
                length(unique(cbind(Biobase::pData(Biobase::phenoData(reactiveBce())),
                                    pred(getStrat(reactiveBce(), input$strategy))
                )[, annot]))
              }))) * 22 - 106)})

      # render the top tab panel
      output$top_tabs <- renderUI({
        tabPanel("Summary", sideBarLayout(
          uiOutput("uiabundance", width = "100%"),
                 plotOutput("pca", width = "100%"))
      )
      })
      
      # render the inside tab panel
      output$inside_tabs <- renderUI({
        tabsetPanel(
          
          tabPanel(
            "Sample distance",
            plotOutput("distance",
                       width = "100%"
            )
          ),
          # tabPanel("Stability", plotOutput("stability", width = "100%")),
          tabPanel("Bicluster members", uiOutput("uiScoreHeatmap"), 
                   plotOutput("score_threshold", 
                              width = "100%")
                   
          ),
          tabPanel("Biomarkers", plotOutput("loadingHeatmap", width = "100%"), 
            plotOutput("plot_biomarkers", width = "100%")),
          id = "main_panel",
          type = "pills"
        )
      })
      
      # plot cluster stability (not implemented yet)
      output$stability <- renderPlot({
        gt <- reactiveStability()
        print(gt)
      })
      reactiveStability <- reactive({
        withProgress(message = "Plotting...", value = 0, {
          set.seed(1234567)
          plotClustStab(reactiveBce(), input$strategy)
        })
      })
      
      # heatmap of scores for all samples
      # using renderUI prevents overlapping plots
      output$uiScoreHeatmap <- renderUI({
        height = reactiveHeatmapHeight300()
        plotOutput("scoreHeatmap", height = height)
        })
      output$scoreHeatmap <- renderPlot({
        gt <- reactiveScoreHeatmap()
        print(gt)
      }, 
      height = function() { reactiveHeatmapHeight300() })
      reactiveScoreHeatmap <- reactive({withProgress(
        message = "Plotting...",
        value = 0, {
          set.seed(1234567)
          phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(reactiveBce())))
          biclustLabels <- intersect(input$annots, names(getStrat(reactiveBce(), input$strategy)))
          heatmapFactor(reactiveBce(), strategy = input$strategy, type = "score",
                         phenoLabels = phenoLabels, 
                         biclustLabels = biclustLabels, ordering = input$sampOrder, 
                         colNames = input$sampNames)
        }
      )})
      
      # Plot of sample scores for one bicluster
      output$score_threshold <- renderPlot({
        reactive_score_threshold()
      })
      reactive_score_threshold <- reactive({
        withProgress(message = "Plotting...", value = 0, {
          set.seed(1234567)
          plotSamples(reactiveBce(), strategy = input$strategy, 
                      bicluster = input$cluster,
                      ordering = input$sampOrder)
        })
      })
      
      # Heatmap of loadings for all features
      output$loadingHeatmap <- renderPlot({
        gt <- reactiveLoadingHeatmap()
        print(gt)
      })
      reactiveLoadingHeatmap <- reactive({withProgress(
        message = "Plotting...",
        value = 0, {
          set.seed(1234567)
          heatmapFactor(reactiveBce(), strategy = input$strategy, type = "loading",
                         ordering = input$featOrder, 
                         colNames = input$featNames)
        }
      )})
      
      # Plot of feature loadings for one bicluster
      output$plot_biomarkers <- renderPlot({
        reactiveMarkers()
      })
      reactiveMarkers <- reactive({
        withProgress(message = "Plotting...", value = 0, {
          set.seed(1234567)
          plotMarkers(reactiveBce(), strategy = input$strategy, 
                      bicluster = input$cluster,
                      ordering = input$featOrder)
        })
      })
      
      # PANEL DESCRIPTIONS
      output$explanation <- renderUI({
        res <- ""
        if (length(input$main_panel) > 0) {
          if ("Abundance" == input$main_panel) {
            res <- HTML(
              "Abundance shows a Z-score heatmap of the original matrix
              of metabolite peaks. Default ordering is distance-based by
              cluster membership. Sidebar options allow reordering
              of rows to highlight samples that are members of the
              currently selected cluster.
              Currently the annotations 'species' and 'tissue' are fictitious."
            )
          }
          # if ("Stability" == input$main_panel) {
          #   res <- HTML(
          #     "Stability index shows how stable each cluster
          #     is accross the selected range of <i>k</i>s.
          #     The stability index varies between 0 and 1, where
          #     1 means that the same cluster appears in every
          #     solution for different <i>k</i>."
          #   )
          # }
          if ("Biomarkers" == input$main_panel) {
            res <- HTML(
              paste0(
                "This is a plot of loading, a statistic that
                quantifies the importance of each feature in distinguishing
                samples that are members of this bicluster."
              )
            )
          }
          if ("Bicluster members" == input$main_panel) {
            res <- HTML(
              "This is a score thresholding plot, showing which
              samples are above the threshold. The user should be
              told which thresholding method was used."
            )
          }
          return(res)
        }
      })
      
      # plotHeightMark <- function() {
      #   return(150 + 10.8 * nrow(values$mark.res))
      # }
      
      # REACTIVE BUTTONS
      # is_biology <- reactive({
      #   return(biology)
      # })
      #
      
      #### REACTIVE DATA #######################################################
      output$selectCluster <- renderUI({
        selectInput("cluster", "Select bicluster:", 
                    choices = names(getStrat(reactiveBce(), input$strategy))
        )
      })
      output$annotPicker <- renderUI({
        classes <- if(inherits(reactiveBce(), "BiclusterExperiment")) {
          c(colnames(Biobase::phenoData(reactiveBce())))
        } else NULL
        selectInput("annots", label = "Annotations",
                    choices = classes,
                    selected = NULL, multiple = TRUE)
      }
      )
      
      # observer for marker genes
      # observe({
      #   if (FALSE) {
      #     # biology) {
      #     # get all marker genes
      #     markers <- organise_marker_genes(object, input$strategy, 
      #                                      as.numeric(input$pValMark), 
      #                                      as.numeric(input$auroc.threshold))
      #     user$n.markers <- nrow(markers)
      #     # get top 10 marker genes of each cluster
      #     markers <- markers_for_heatmap(markers)
      #     clusts <- unique(markers[, 1])
      #     if (is.null(clusts)) {
      #       clusts <- "None"
      #     }
      #     user$mark.res <- markers
      #     updateSelectInput(session, "cluster", choices = clusts)
      #   } else {
      #     user$n.markers <- 0
      #   }
      # })
      # 
      # 
      # output$has_biomarkers <- reactive({
      #   return(TRUE)
      # })
      
      # stop App on closing the browser
      session$onSessionEnded(function() {
        stopApp()
      })
      # 
      # outputOptions(output, "has_biomarkers", suspendWhenHidden = FALSE)
    },
    
    # launch App in a browser
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})