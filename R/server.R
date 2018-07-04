function(input, output, session) {
  dep <- reactiveValues(
    pheatmap = requireNamespace("pheatmap", quietly = TRUE),
    plaid = requireNamespace("biclust", quietly = TRUE),
    spectral = requireNamespace("biclust", quietly = TRUE)
  )
  
  values <- reactiveValues(
    # Initialize as the user's BiclusterExperiment.
    bce = if(is.null(userBce)) {
      NULL
    } else userBce,
    rawmat = if(!is.null(userBce)) as.matrix(userBce) else NULL,
    strategy = NULL
  )
  reactiveSeparator <- observeEvent(input$input_df, {
    # If the file input changes, try automatically changing separator to ","
    if(substr(input$input_df$name, nchar(input$input_df$name) - 2,
              nchar(input$input_df$name)) == "csv") {
      updateTextInput(session, "sepchar", value = ",")
    }
  }
  )
  
  observeEvent({ # read the raw input as a data.frame
    input$input_df
    input$row_names
    input$header
    input$sepchar
    input$quotechar
    input$decchar
    input$skiplines
    input$sampleCols
  }, 
  { # If the user has not chosen a file yet, look for a BiclusterExperiment in the execution environment of biclusterGUI
    if(is.null(input$input_df)) {
      shinyjs::disable("row_names")
      shinyjs::disable("header")
      shinyjs::disable("skiplines")
      shinyjs::disable("sepchar")
      shinyjs::disable("quotechar")
      shinyjs::disable("decchar")
      validate(need(exists("user") && !is.null(userBce), 
                    "You may import your dataset."))
      as.data.frame(t(as.matrix(userBce)))
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
      withProgress({
        incProgress(1/8)
        rawmat <- read.table(input$input_df$datapath, header = input$header,
                             row.names = rows, fill = TRUE, comment.char = "",
                             sep = input$sepchar, quote = input$quotechar, dec = input$decchar, skip = input$skiplines)
        incProgress(6/8)
        rawmat <- as.matrix(rawmat)
        if(input$sampleCols) rawmat <- t(rawdf)
        incProgress(1/8)
      }, message = "Parsing data...", value = 0)
      values$rawmat <- rawmat
    }
  })
  
  # Always keep the BiclusterExperiment synchronized with the user's fi
  observeEvent(values$rawmat, {
    values$bce <- BiclusterExperiment(values$rawmat)
  })
  
  output$dt <- DT::renderDT({
    return(values$rawmat)
  }, server = FALSE, selection = "none")
  reactiveDF <- reactive({
    values$rawmat
  })
  
  # helper functions
  reactiveHeatmapHeight500 <- reactive({ 
    max(500, length(input$annots) * 33 +
          sum(unlist(lapply(input$annots, function(annot) {
            length(unique(cbind(Biobase::pData(Biobase::phenoData(values$bce)),
                                pred(getStrat(values$bce, input$strategy))
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
  
  output$annotPicker <- renderUI({
    classes <- if(inherits(values$bce, "BiclusterExperiment")) {
      c(colnames(Biobase::phenoData(values$bce)))
    } else NULL
    selectInput("annots", label = "Annotations",
                choices = classes,
                selected = NULL, multiple = TRUE)
  }
  )
  
  # plot abundance heatmap (original data)
  output$uiabundance <- renderUI({
    validate(need(inherits(values$bce, "BiclusterExperiment"), "You may import your dataset."))
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
    bce <- values$bce
    
    withProgress(message = "Plotting...", value = 0, {
      set.seed(1234567)
      if(input$logBase == "e") {
        logBase <- exp(1)
      } else {
        logBase <- as.numeric(input$logBase)
      }
      phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(bce)))
      # biclustLabels <- if(length(names(bce)) > 0) {
      #   intersect(input$annots, names(getStrat(bce, input$strategy)))
      plot(bce, logBase = logBase, phenoLabels = phenoLabels,
           ordering = if(input$heatmapReorder) "distance" else "input",
           strategy = input$strategy,
           rowNames = input$sampNames, colNames = input$featNames)
    })
  })
  
  # Plot of samples along first two PCs
  output$pca <- renderPlot({ 
    validate(need(inherits(values$bce, "BiclusterExperiment"), "You may import your dataset."))
    pca(values$bce)
  }
  )
  
  # Euclidean Distance between samples (applied to raw data...might be good to scale first?)
  output$distance <- renderPlot({
    validate(need(inherits(values$bce, "BiclusterExperiment"), "You may import your dataset."))
    gt <- reactiveDistance()
    print(gt)
  }, height = function() {reactiveHeatmapHeight500()})
  reactiveDistance <- reactive({
    validate(need(
      all(input$annot.rows %in% colnames(phenoData(values$bce))),
      paste0(
        "\nPlease choose annotations from: ",
        paste(colnames(phenoData(values$bce)), collapse = ", ")
      )
    ))
    bce <- values$bce
    withProgress(message = "Plotting...", value = 0, {
      set.seed(1234567)
      phenoLabels <- intersect(input$annots, colnames(Biobase::phenoData(bce)))
      # biclustLabels <- intersect(input$annots, names(getStrat(bce, input$strategy)))
      biclustLabels <- NULL
      distType <- if(input$distType == "Euclidean distance") {
        "euclidean"} else {"pearson"}
      plotDist(bce, distType = distType, phenoLabels = phenoLabels, 
               biclustLabels = biclustLabels, 
               ordering = if(input$distReorder) "distance" else "input", 
               strategy = input$strategy, 
               rowColNames = input$sampNames
      )
    })
  })
  
  #### BICLUSTER ####
  
  output$kSlider <- renderUI({
    bce <- values$bce
    withProgress({
      maxk <- if(is.null(bce)) { 2 }
      else {
        m <- as.matrix(bce)
        min(nrow(m), ncol(m))
      }
      sliderInput("k", "Number biclusters", value = 2, min = 1,
                  max = maxk, step = 1)
    },
    message = "Loading...",
    value = 0)
  }
  )
  
  observeEvent(input$bicluster, {
    validate(need(inherits(values$bce, "BiclusterExperiment"), "You may import your dataset."))
    bce <- values$bce
    # Look for an existing BiclusterStrategy with the same parameters
    stratName <- name(list(bca = capitalize(input$algo), 
                           sta = capitalize("otsu"), lta = capitalize("otsu"), k = input$k))
    matchStrats <- which(names(bce) == stratName)
    if(length(matchStrats) == 0) {
      withProgress({
        bce <- addStrat(bce, k = input$k, method = tolower(input$algo))
        newStrat <- names(bce)[length(names(bce))]
        algo <- strsplit(newStrat, split = " | ")[[1]]
        if(algo != input$algo) {
          updateSelectInput(session, inputId = "algo", selected = algo)
          showNotification(
            paste(input$algo, "failed on your dataset, so the", algo,
                  "algorithm was used instead."), duration = 5)
        }
      }, message = "Biclustering...", value = 0.2)
    } else {
      values$strategy <- getStrat(bce, matchStrats[1])
    }
    shinyjs::disable("bicluster")
  })
  observeEvent({values$bce
    input$algo
    input$k
  }, shinyjs::enable("bicluster")
  )
  
  #### Scores ####
  output$selectCluster <- renderUI({
    stratName <- name(list(bca = capitalize(input$algo), 
                           sta = capitalize("otsu"), 
                           lta = capitalize("otsu"), k = input$k))
    selectInput("cluster", "Select bicluster:", 
                choices = names(getStrat(values$bce, stratName))
    )
  })
  
  reactiveHeatmapHeight300 <- reactive({
    bce <- values$bce
    max(300, length(input$annots) * 33 +
          sum(unlist(lapply(input$annots, function(annot) {
            length(unique(cbind(Biobase::pData(Biobase::phenoData(bce)),
                                pred(getStrat(bce, input$strategy))
            )[, annot]))
          }))) * 22 - 106)})
  
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
  reactiveScoreHeatmap <- reactive({
    bce <- values$bce
    withProgress(
      message = "Plotting...",
      value = 0, {
        set.seed(1234567)
        phenoLabels <- intersect(input$annots, 
                                 colnames(Biobase::phenoData(bce)))
        biclustLabels <- intersect(input$annots, 
                                   names(getStrat(values$bce, 
                                                  input$strategy)))
        heatmapFactor(bce, strategy = input$strategy, type = "score",
                      phenoLabels = phenoLabels, 
                      biclustLabels = biclustLabels, 
                      ordering = if(input$scoreReorder) { "cluster" }
                      else { "input" }, 
                      colNames = input$sampNames)
      })
  })
  
  # Plot of sample scores for one bicluster
  output$score_threshold <- renderPlot({
    reactive_score_threshold()
  })
  reactive_score_threshold <- reactive({
    withProgress(message = "Plotting...", value = 0, {
      set.seed(1234567)
      plotSamples(values$bce, strategy = input$strategy, 
                  bicluster = input$cluster,
                  ordering = if(input$sampOrder) { "input" }
                  else { "distance" })
    })
  })
  
  
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
      plotClustStab(values$bce, input$strategy)
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
      heatmapFactor(values$bce, strategy = input$strategy, type = "loading",
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
      plotMarkers(values$bce, strategy = input$strategy, 
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
    }