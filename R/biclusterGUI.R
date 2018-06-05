#' @include BiclusterExperiment.R
NULL

#' Explore a BiclusterExperiment
#'
#' Opens a GUI to visualize analysis results contained in a BiclusterExperiment
#' object
#'
#' @export
setGeneric("biclusterGUI", signature = "obj", function(obj) {
  standardGeneric("biclusterGUI")
})

# clusters must be a named list of matrices
#' @describeIn biclusterGUI Open GUI for a BiclusterExperiment
#' @importFrom shiny HTML actionButton animationOptions checkboxInput checkboxGroupInput column div downloadHandler downloadLink eventReactive fluidPage fluidRow h4 headerPanel htmlOutput need observe observeEvent p plotOutput reactiveValues renderPlot renderUI selectInput shinyApp sliderInput stopApp tabPanel tabsetPanel uiOutput updateSelectInput validate wellPanel withProgress conditionalPanel reactive outputOptions tags radioButtons downloadButton sidebarLayout sidebarPanel mainPanel
setMethod("biclusterGUI", c(obj = "BiclusterExperiment"), function(obj) {
  ## define UI parameters
  plotHeight <- 600
  plotHeightSmall <- 300
  
  ## define server global variables
  values <- reactiveValues()
  
  shinyApp(
    ui = fluidPage(
      #### UI ##################################################################
      tags$head(tags$style(
        HTML(".shiny-output-error-validation {
             color: red;
}")
      )),
      fluidRow(column(
        12,
        HTML("<br>")
      )),
      
      sidebarLayout(
        # User Control panel
        sidebarPanel(
          tags$head(tags$style("#loadingHeatmap{height:30vh !important;}")),
          wellPanel(
            # "Select data" panel. Enables user to select the BiclusterStrategy
            # and which bicluster in that strategy.
            # TODO add file I/O buttons.
            h4("Select data"),
            selectInput(inputId = "strategy", 
                        label = "Biclustering strategy",
                        choices = names(obj))
          ),
            # FIXME: Try changing the default to the BCV cluster number?
          wellPanel(fluidRow(
            column(
              12,
              HTML("<font size = 2>"),
              h4("Format"),
              # Log transform the original data (only for
              # abundance. We expect the user to be
              # responsible for homoscedascity,
              # normalization, etc.
              conditionalPanel("input.main_panel == 'Abundance'",
                               selectInput("logBase", "Log-transform",
                                           choices = list("None" = 0, 
                                                          "log[2]" = 2, 
                                                          "ln" = "e", 
                                                          "log[10]" = 10))),
              # Formula for input data
              conditionalPanel("input.main_panel == 'Sample distance'",
                               selectInput("distType", "Formula",
                                           choices = c("Euclidean distance",
                                                       "1 - Pearson"))),
              # Choose bicluster to highlight (for score and loading plots)
              conditionalPanel("input.main_panel == 'Bicluster members' ||
                                input.main_panel == 'Biomarkers'",
                               uiOutput("selectCluster")),
              # Choose sample ordering
              conditionalPanel("input.main_panel != 'Biomarkers'",
                               selectInput(
                                 "sampOrder",
                                 "Sample ordering",
                                 choices = c(
                                   "Same as input" = "input",
                                   "Euclidean dist by abundance" = "distance",
                                   "Euclidean dist by bicluster membership" = "cluster"
                                 )
                               )),
              # Choose feature ordering (only for Biomarkers panel)
              conditionalPanel("input.main_panel == 'Biomarkers'",
                               selectInput(
                                 "featOrder",
                                 "Feature ordering",
                                 choices = c(
                                   "Same as input" = "input"
                                   # ,
                                   # "Euclidean dist by bicluster membership" = "cluster"
                                 )
                               )),
              # Pick annotations to show along plot axes
              conditionalPanel("input.main_panel != 'Stability' && 
                               input.main_panel != 'Biomarkers'",
                               uiOutput("annotPicker")),
              # Show sample names?
              conditionalPanel("input.main_panel != 'Stability' &&
                                 input.main_panel != 'Biomarkers'",
                               checkboxInput("sampNames", "Show sample names")),
              # Show feature names?
              conditionalPanel("input.main_panel == 'Abundance' ||
input.main_panel == 'Biomarkers'",
                               checkboxInput("featNames", "Show feature names")),
              HTML("</font>")
            )
          )),
          width = 3
        ),
        
        # Data panel and description column
        mainPanel(
          column(
            9,
            uiOutput("mytabs")
          ),
          column(
            3,
            h4("Panel description"),
            htmlOutput("explanation")
          ),
          width = 9
        )
      )
    ),
    
    #### SERVER ##################################################################
    server = function(input, output, session) {
      # render the top tab panel
      output$mytabs <- renderUI({
        tabsetPanel(
          tabPanel("Abundance", plotOutput("abundance", width = "100%")),
          tabPanel(
            "Sample distance",
            plotOutput("distance",
                       width = "100%"
            )
          ),
          tabPanel("Stability", plotOutput("stability", width = "100%")),
          tabPanel("Bicluster members", plotOutput("scoreHeatmap", 
                                                   width = "100%"), 
                   plotOutput("score_threshold", 
                              width = "100%")
                   
          ),
          tabPanel("Biomarkers", plotOutput("loadingHeatmap", width = "100%"), 
            plotOutput("plot_biomarkers", width = "100%")),
          id = "main_panel",
          type = "pills"
        )
      })
      #### REACTIVE PANELS #######################################################
      # plot abundance heatmap (original data)
      reactiveHeatmapHeight500 <- reactive({ max(500, length(input$annots) * 33 +
          sum(unlist(lapply(input$annots, function(annot) {
            length(unique(cbind(Biobase::pData(Biobase::phenoData(obj)),
                                pred(getStrat(obj, input$strategy))
            )[, annot]))
          }))) * 22 - 106)})
      
      reactiveHeatmapHeight300 <- reactive({ max(300, length(input$annots) * 33 +
                                                   sum(unlist(lapply(input$annots, function(annot) {
                                                     length(unique(cbind(Biobase::pData(Biobase::phenoData(obj)),
                                                                         pred(getStrat(obj, input$strategy))
                                                     )[, annot]))
                                                   }))) * 22 - 106)})
      output$abundance <- renderPlot({
        reactive_abundance()
      }, height = function() {reactiveHeatmapHeight500()})
      
      
      reactive_abundance <- reactive({
        withProgress(message = "Plotting...", value = 0, {
          set.seed(1234567)
          if(input$logBase == "e") {
            logBase <- exp(1)
          } else {
            logBase <- as.numeric(input$logBase)
          }
          phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(obj)))
          biclustLabels <- intersect(input$annots, names(getStrat(obj, input$strategy)))
          plot(obj, logBase = logBase, phenoLabels = phenoLabels, biclustLabels = biclustLabels, ordering = input$sampOrder, strategy = input$strategy, 
               rowNames = input$sampNames, colNames = input$featNames)
        })
      })
      
      # Euclidean Distance between samples (applied to raw data...might be good to scale first?)
      output$distance <- renderPlot({
        reactiveDistance()
      }, height = function() {reactiveHeatmapHeight500()})
      
      reactiveDistance <- reactive({
        validate(need(
          all(input$annot.rows %in% colnames(phenoData(obj))),
          paste0(
            "\nPlease choose annotations from: ",
            paste(colnames(phenoData(obj)), collapse = ", ")
          )
        ))
        withProgress(message = "Plotting...", value = 0, {
          set.seed(1234567)
          phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(obj)))
          biclustLabels <- intersect(input$annots, names(getStrat(obj, input$strategy)))
          distType <- if(input$distType == "Euclidean distance") {
            "euclidean"} else {"pearson"}
          plotDist(obj, distType = distType, phenoLabels = phenoLabels, 
                   biclustLabels = biclustLabels, 
                   ordering = input$sampOrder, strategy = input$strategy, 
                   rowColNames = input$sampNames
          )
        })
      })
      
      # plot cluster stability (not implemented yet)
      output$stability <- renderPlot({
        reactiveStability()
      })
      
      reactiveStability <- reactive({
        withProgress(message = "Plotting...", value = 0, {
          set.seed(1234567)
          plotClustStab(obj, input$strategy)
        })
      })
      
      # heatmap of scores for all samples
      output$scoreHeatmap <- renderPlot({
        reactiveScoreHeatmap()
      }, height = function() {reactiveHeatmapHeight300()})
      reactiveScoreHeatmap <- reactive({withProgress(
        message = "Plotting...",
        value = 0, {
          set.seed(1234567)
          phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(obj)))
          biclustLabels <- intersect(input$annots, names(getStrat(obj, input$strategy)))
          heatmapFactor(obj, strategy = input$strategy, type = "score",
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
          plotSamples(obj, strategy = input$strategy, 
                      bicluster = input$cluster,
                      ordering = input$sampOrder)
        })
      })
      
      # Heatmap of loadings for all features
      output$loadingHeatmap <- renderPlot({
        reactiveLoadingHeatmap()
      })
      reactiveLoadingHeatmap <- reactive({withProgress(
        message = "Plotting...",
        value = 0, {
          set.seed(1234567)
          heatmapFactor(obj, strategy = input$strategy, type = "loading",
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
          plotMarkers(obj, strategy = input$strategy, 
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
              "ADD GUI BUTTONS FOR IMPORTING DATA, SWITCHING TO PCA, BIMAX, PLAID etc. 
              Cause annotation colors to not get shuffled. Abundance shows a Z-score heatmap of the original matrix
              of metabolite peaks. Default ordering is distance-based by
              cluster membership. Sidebar options allow reordering
              of rows to highlight samples that are members of the
              currently selected cluster.
              Currently the annotations 'species' and 'tissue' are fictitious."
            )
          }
          if ("Stability" == input$main_panel) {
            res <- HTML(
              "Stability index shows how stable each cluster
              is accross the selected range of <i>k</i>s.
              The stability index varies between 0 and 1, where
              1 means that the same cluster appears in every
              solution for different <i>k</i>."
            )
          }
          if ("Biomarkers" == input$main_panel) {
            res <- HTML(
              paste0(
                "ADD ACCURACY, RECOVERY, etc. on
                PLOT. ALLOW TO HIGHLIGHT KNOWN FEATURES. This is a plot of loading, a statistic that
                quantifies the importance of each feature in distinguishing
                samples that are members of this bicluster.
                Biomarkers are contained in the generated msNMF object."
              )
            )
          }
          if ("Bicluster members" == input$main_panel) {
            res <- HTML(
              "ADD ACCURACY, RECOVERY, etc. on
              PLOT. ALLOW TO HIGHLIGHT KNOWN SAMPLE GROUPS. This is a score thresholding plot, showing which
              samples are above the threshold. The user should be
              told which thresholding method was used.
              Currently the values on the y-axis are fictitious
              score values for fictitious bicluster 1."
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
      
      #### REACTIVE INPUTS #######################################################
      output$selectCluster <- renderUI({
        selectInput("cluster", "Select bicluster:", 
                    choices = names(getStrat(obj, input$strategy))
        )
      })
      output$annotPicker <- renderUI({
        selectInput("annots", label = "Annotations",
                    choices = c(colnames(Biobase::phenoData(obj)), 
                                names(getStrat(obj, input$strategy))),
                    selected = NULL, multiple = TRUE)
      }
      )
      
      # observer for marker genes
      observe({
        if (FALSE) {
          # biology) {
          # get all marker genes
          markers <- organise_marker_genes(object, input$strategy, 
                                           as.numeric(input$pValMark), 
                                           as.numeric(input$auroc.threshold))
          values$n.markers <- nrow(markers)
          # get top 10 marker genes of each cluster
          markers <- markers_for_heatmap(markers)
          clusts <- unique(markers[, 1])
          if (is.null(clusts)) {
            clusts <- "None"
          }
          values$mark.res <- markers
          updateSelectInput(session, "cluster", choices = clusts)
        } else {
          values$n.markers <- 0
        }
      })
      
      
      output$has_biomarkers <- reactive({
        return(TRUE)
      })
      
      # stop App on closing the browser
      session$onSessionEnded(function() {
        stopApp()
      })
      
      outputOptions(output, "has_biomarkers", suspendWhenHidden = FALSE)
    },
    
    # launch App in a browser
    options = list(launch.browser = TRUE, fullstacktrace = TRUE)
  )
})