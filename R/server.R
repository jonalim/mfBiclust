function(input, output, session) {
  
  params <- list(
    biclusterargs = if(debug) {
      list(maxIter = 100L, shuffle = 1, withinVar = 10)
      } else list(),
    bcvMaxIter = if(debug) 3L else 100L, # how many iterations of bcv to run
    annotateBiclusters = if(debug) 3L else NA # how many biclusters to annotate
  )

  #### REACTIVE VALUES #### 
  dep <- reactiveValues(
    pheatmap = requireNamespace("pheatmap", quietly = TRUE),
    plaid = requireNamespace("biclust", quietly = TRUE),
    spectral = requireNamespace("biclust", quietly = TRUE)
  )
  
  values <- reactiveValues(
    # Initialize as the user's BiclusterExperiment.
    bce = if(inherits(userBce, "BiclusterExperiment")) userBce else NULL,
    strategy = if(inherits(userBce, "BiclusterExperiment")) {
      if(length(strategies(userBce)) > 0) {
        getStrat(userBce, 1)
      }
    } else NULL,
    zoom = c(0, 1, 0, 1),
    bcvRes = numeric(0),
    bcvBest = 0,
    bcvValid = FALSE,
    # a named list containing results for each bicluster
    goRes = numeric(0),
    goValid = FALSE
  )
  
  #### DATA I/O ####
  reactiveSeparator <- observeEvent(input$input_df, {
    # If the file input changes, try automatically changing separator to ","
    if(substr(input$input_df$name, nchar(input$input_df$name) - 2,
              nchar(input$input_df$name)) == "csv") {
      updateTextInput(session, "sepchar", value = ",")
    }
  }
  )
  
  rawmat <- reactive({ # If the user has not chosen a file yet, look for a BiclusterExperiment in the execution environment of biclusterGUI
    if(is.null(input$input_df)) {
      shinyjs::disable("row_names")
      shinyjs::disable("header")
      shinyjs::disable("skiplines")
      shinyjs::disable("sepchar")
      shinyjs::disable("quotechar")
      shinyjs::disable("decchar")
      validate(need(inherits(userBce, "BiclusterExperiment"), 
                    "You may import your dataset."))
      t(as.matrix(userBce))
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
      rawmat
    }
  })
  # 
  # customRowNames <- observeEvent(
  #   input$customRowNames,
  #   {
  #     if(nchar(input$customRowNames) > 0) {
  #       rn <- scan(text = input$customRowNames, what = "",
  #            quiet = TRUE)
  #       if(inherits(rawmat, "matrix")) {
  #         try(rownames(rawmat) <- rn)
  #       }
  #       if(inherits(values$bce, BiclusterStrategy)) {
  #         
  #              } else return(character(0))
  #              }
  # )
  
  # When the raw data matrix changes, trigger other reactive updates
  observeEvent(rawmat(), {
    values$bce <- BiclusterExperiment(rawmat())
    bcvValid <- FALSE
  })
  
  # Render an interactive data.table for the user's benefit
  output$dt <- DT::renderDT({
    return(rawmat())
  }, server = FALSE, selection = "none")
  
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
  
  #### BICLUSTER TAB ####
  # FIXME need to initialize the sliderInput before even navigating to this tab.
  # That way, BCV and Bicluster can modify the sliderInput.
  # TODO initialize the slidreInput in the UI, then use an observer on values$bce
  # to modify the max and value ASAP
  output$kSlider <- renderUI({
    bce <- values$bce
    val <- input$k
    if(is.null(val)) {
      # Initialize everything on the Bicluster tab
      if(inherits(values$strategy, "BiclusterStrategy")) {
        # Assume a BiclusterStrategy has a method listed in biclusterUI.R
        updateSelectInput(session, "algo", selected = method(values$strategy))
        val <- nclust(values$strategy)
      } else val <- 2
    }
    maxk <- if(is.null(bce)) { 2 } else {
      m <- as.matrix(bce)
      min(nrow(m), ncol(m))
    }
    sliderInput("k", "Number biclusters", min = 1, 
                value = val,
                max = maxk, step = 1)
  }
  )
  
  # observeEvent({values$bce
  #   input$k
  # },
  # {if(!(inherits(values$bce, "BiclusterExperiment") && inherits(input$k, "numeric"))) {
  #   shinyjs::disable("bicluster")
  # } else {
  #   shinyjs::enable("bicluster")
  # }
  # }
  # )
  
  observeEvent(input$bicluster, {
    validate(need(inherits(values$bce, "BiclusterExperiment"), "You may import your dataset."))
    validate(need(inherits(input$k, "integer"), ""))
    shinyjs::disable("bicluster")
    try(removeNotification("overlapNotif"))
    
    bce <- values$bce
    # Look for an existing BiclusterStrategy with the same parameters
    stratName <- name(list(bca = capitalize(input$algo), 
                           sta = capitalize("otsu"), lta = capitalize("otsu"), k = input$k))
    matchStrats <- which(names(bce) == stratName)
    if(length(matchStrats) == 0) {
      withProgress({
        if("withinVar" %in% names(params$biclusterargs)) {
          params$biclusterargs$withinVar <- params$biclusterargs$withinVar * 
            nrow(bce)
        }
        # append any optional debug-mode arguments
        withCallingHandlers(
          bce <- do.call(addStrat,
                         c(bce = bce, k = input$k, method = tolower(input$algo),
                                     duplicable = TRUE, params$biclusterargs)),
        warning = function(w) {
          # In case less than the requested number of biclusters was found
          showNotification(w$message, duration = NULL)
        })
        newStrat <- names(bce)[length(names(bce))]
        algo <- strsplit(newStrat, split = " | ")[[1]][1]
        if(algo != input$algo) {
          updateSelectInput(session, inputId = "algo", selected = algo)
          showNotification(
            paste(input$algo, "failed on your dataset, so the", algo,
                  "algorithm was used instead."), duration = 5)
        }
        # Whatever the new BiclusterStrategy is, set it as the active strategy
        values$strategy <- getStrat(bce, newStrat)
        values$bce <- bce
      }, message = "Biclustering...", value = 0.2)
    } else {
      values$strategy <- getStrat(bce, matchStrats[1])
    }
  })
  
  # Manage everything that happens when the parameters are changed
  observeEvent({values$bce
    input$algo
    input$k
  }, {
    # Disable the Run button if the data isn't in the right format
    if(!(inherits(values$bce, "BiclusterExperiment") && inherits(input$k, "integer"))) {
      shinyjs::disable("bicluster") 
    } else {
      # Enable the run button if the Bicluster strategy hasn't been run yet.
      bce <- values$bce
      stratName <- name(list(bca = capitalize(input$algo), 
                             sta = capitalize("otsu"), lta = capitalize("otsu"), k = input$k))
      matchStrats <- which(names(bce) == stratName)
      if(length(matchStrats) == 0) {
        shinyjs::enable("bicluster")
      } else {
        # If the BiclusterStrategy has been run already, activate it and disable
        # the Run button
        shinyjs::disable("bicluster")
        values$strategy <- strategies(bce)[[stratName]]
      }
    }
  }
  )
  
  #### Summary ####
  # Reset zoom when the user's input data changes or upon double-click  
  observeEvent({rawmat()
    input$image_dblclick
  }, {
    mat <- rawmat()
    
    if(inherits(mat, "matrix") && mode(mat) == "numeric") {
      values$zoom <- c(0, 1, 0, 1)
    }
  })
  
  # Re-render the bicluster plot when raw data OR BiclusterStrategy OR selected biclusters changes
  imageArr <- reactive({
    # render the whole heatmap
    mat <- rawmat()
    validate(need(inherits(mat, "matrix") && mode(mat) == "numeric", ""))
    m <- mat / max(mat)
    width <- ncol(m)
    height <- nrow(m)
    # Convert the vector to an array with 3 planes for R, G, B
    arr <- array(c(m, m, m), dim = c(height, width, 3))
    
    # Add bicluster outlines 
    bcs <- values$strategy
    if(!is.null(bcs)) {
      cols <- hcl(h = seq(0, (nclust(bcs) - 1) / (nclust(bcs)),
                          length = nclust(bcs)) * 360, c = 100, l = 65,
                  fixup = TRUE)
      cols <- col2rgb(cols) # RGB as rows; each column a different color
      # FIXME: Allow to select a subset of biclusters. Allow bicluster ID when
      # hovering mouse over.
      lapply(seq_len(nclust(bcs)), function(bicluster) {
        yrange <- which(loading(bcs)[bicluster, ] > bcs@loadingThresh[1, 1])
        xrange <- which(pred(bcs)[, bicluster])
        arr[xrange, yrange, 1] <<- cols["red", bicluster] / 255
        arr[xrange, yrange, 2] <<- cols["green", bicluster] / 255
        arr[xrange, yrange, 3] <<- cols["blue", bicluster] / 255
      })
    }
    arr
  })
  
  output$image1 <- renderImage({
    arr <- imageArr()
    validate(need(inherits(arr, "array") && 
                    mode(arr) == "numeric" &&
                    !is.null(dim(arr)),
                  "Data for heatmap has not been parsed yet. Try switching revisiting the 'Data' tab momentarily."))
    temp <- tempfile(fileext = ".png")
    # if(inherits(arr, "array") && mode(arr) == "numeric") {
    
    # subset the PNG and re-render
    rangeY <- c(floor(values$zoom[3] * dim(arr)[1]) + 1, round(values$zoom[4] * dim(arr)[1]))
    rangeX <- c(floor(values$zoom[1] * dim(arr)[2]) + 1, round(values$zoom[2] * dim(arr)[2]))
    # bug; it's possible to zoom in so much that rangeY[2] - rangeY[1]< 1 similarly for rangeX
    
    arr <- arr[rangeY[1]:rangeY[2], rangeX[1]:rangeX[2], ]
    # } else {
    #   arr <- array(1, dim = c(1, 1, 1)) # one white pixel
    #   
    # }
    tryCatch(png::writePNG(arr,
                           target = temp), error = function(e) {
                             browser()}
    )
    # Return a list containing information about the zoomed image
    list(
      src = temp,
      contentType = "image/png",
      title = "Click and drag to zoom in; double-click to reset"
    )
  })
  
  # If the user uses the brush, zoom in. Clearing the brush simply allows the
  # GUI to await another brush input.
  observeEvent(input$image_brush, {
    brush <- input$image_brush
    zoom <- values$zoom
    if(!is.null(brush)) {
      # contains fractions of the whole image
      values$zoom <- c(zoom[1] + (zoom[2] - zoom[1]) * brush$xmin / session$clientData$output_image1_width, 
                       zoom[1] + (zoom[2] - zoom[1]) * brush$xmax / session$clientData$output_image1_width,
                       zoom[3] + (zoom[4] - zoom[3]) * brush$ymin / session$clientData$output_image1_height,
                       zoom[3] + (zoom[4] - zoom[3]) * brush$ymax / session$clientData$output_image1_height)
      shinyjs::runjs("document.getElementById('image1_brush').remove()")
    }
  })
  
  # Warn if biclusters overlap
  overlapWarn <- observeEvent(values$strategy,
                              {
                                validate(need(inherits(values$strategy, "BiclusterStrategy"), ""))
                                bcs <- values$strategy
                                bc <- threshold(loading(bcs), MARGIN = 1,
                                                loadingThresh(bcs))
                                # Create lists of rows and columns contained in biclusters
                                rc <- biclusterMatrix2List(pred(values$strategy), bc)
                                biclusterRows <- rc[[1]]
                                biclusterCols <- rc[[2]]
                                
                                # Check for overlaps with any other biclusters. Since a bicluster
                                # always overlaps with itself completely, test sum(l) == 1
                                overlaps <- overlap(biclusterRows, biclusterCols)
                                nonOverlap <- sapply(overlaps, function(l) sum(l == 1))
                                if(!all(nonOverlap)) {
                                  showNotification(id = "overlapNotif", paste(
                                    do.call(paste, c(as.list(names(bcs)[!nonOverlap]), sep = ", ")),
                                    "overlap. Please interpret the heatmap annotations with care."), 
                                    duration = NULL)
                                }
                              })
  
  
  #### Scores ####
  output$scorePanel <- renderUI({
    fluidRow(
      column(9,
             plotOutput("scoreHeatmap", width = "100%"), 
             plotOutput("scorePlot", width = "100%")),
      column(3,
             uiOutput("scoreBicluster"),
             checkboxInput("scoreReorder", "Reorder"),
             checkboxInput("sampNames", "Sample names"))
    )
  })
  
  output$scoreBicluster <- renderUI({
    choices <- if(inherits(values$strategy, "BiclusterStrategy")) {
      names(values$strategy)
    } else { list() }
    selectInput("scoreBicluster", "Select bicluster:",
                choices = choices)
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
  output$scoreHeatmap <- renderUI({
    height = reactiveHeatmapHeight300()
    plotOutput("scoreHeatmap", height = height)
  })
  output$scoreHeatmap <- renderPlot({
    gt <- reactiveScoreHeatmap()
    print(gt)
  },
  height = function() { reactiveHeatmapHeight300() })
  reactiveScoreHeatmap <- reactive({
    validate(need(inherits(values$strategy, "BiclusterStrategy"), "Please run biclustering first."))
    bce <- values$bce
    withProgress(
      message = "Plotting...",
      value = 0, {
        set.seed(1234567)
        # Create the annotation selector
        phenoLabels <- intersect(input$annots,
                                 colnames(Biobase::phenoData(bce)))
        biclustLabels <- intersect(input$annots,
                                   names(values$strategy))
        heatmapFactor(bce = bce, bcs = values$strategy, type = "score",
                      phenoLabels = phenoLabels,
                      biclustLabels = biclustLabels,
                      ordering = if(input$scoreReorder) { "cluster" }
                      else { "input" },
                      colNames = input$sampNames)
      })
  })
  
  # Plot of sample scores for one bicluster
  output$scorePlot <- renderPlot({
    scorePlotHelper()
  })
  scorePlotHelper <- reactive({
    validate(need(inherits(values$strategy, "BiclusterStrategy"), ""))
    withProgress(message = "Plotting...", value = 0, {
      set.seed(1234567)
      plotSamples(values$bce, strategy = values$strategy,
                  bicluster = input$scoreBicluster,
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
  
  #### Loading #### 
  
  #FIXME Invalidate these panels when user changes kSlider, but hasn't clicked
  #"Run" yet
  # output$loadingPanel <- renderUI({
  #   
  # })
  # 
  output$loadingBicluster <- renderUI({
    choices <- if(inherits(values$strategy, "BiclusterStrategy")) {
      names(values$strategy)
    } else { list() }
    selectInput("loadingBicluster", "Select bicluster:",
                choices = choices)
  })
  
  # Heatmap of loadings for all features
  output$loadingHeatmap <- renderPlot({
    gt <- reactiveLoadingHeatmap()
    print(gt)
  })
  reactiveLoadingHeatmap <- reactive({
    validate(need(inherits(values$strategy, "BiclusterStrategy"), "Please bicluster first."))
    withProgress(
      message = "Plotting...",
      value = 0, {
        set.seed(1234567) # FIXME: use duplicable()
        heatmapFactor(bce = values$bce, bcs = values$strategy, type = "loading",
                      ordering = input$featOrder,
                      colNames = input$featNames)
      }
    )})
  
  # Plot of feature loadings for one bicluster
  output$plot_biomarkers <- renderPlot({
    validate(need(inherits(values$strategy, "BiclusterStrategy") &&
                    !is.null(input$loadingBicluster), ""))
    reactiveMarkers()
  })
  reactiveMarkers <- reactive({
    withProgress(message = "Plotting...", value = 0, {
      set.seed(1234567)
      plotMarkers(values$bce, strategy = name(values$strategy),
                  bicluster = input$loadingBicluster,
                  ordering = "input")
    })
  })
  
  # Get gene list for selected bicluster
  output$biclusterGeneList <- renderText({
    bce <- values$bce
    bcs <-  values$strategy
    bicluster <- input$loadingBicluster
    validate(need(inherits(bce, "BiclusterExperiment") && 
                    !is.null(input$loadingBicluster) &&
                    inherits(bcs, "BiclusterStrategy"), ""))

    geneI <- which(loading(bcs)[bicluster, ] > bcs@loadingThresh[bicluster, 1])
    genes <- featureNames(bce)[geneI]
    paste(unlist(genes), collapse = "\n")
  })
  
  output$biclusterGeneListLabel <- renderUI({
    bicluster <- input$loadingBicluster
    desc <- if(is.null(bicluster)) { "Markers:" } else {
      paste(bicluster, "markers:")
    }
    tags$h6(desc)
  })
  
  #### BI-CROSS-VALIDATION ####
  guiBCV <- function() {
    if(values$bcvValid == FALSE) {
      withProgress(
        message = "Performing bi-cross-validation...",
        value = 0, max = params$bcvMaxIter, {
          res <- withCallingHandlers({
            shinyjs::html("bcvtext", "")
            suppressWarnings(
              auto_bcv(Y = rawmat(), ks = seq_len(nrow(rawmat())), 
                       bestOnly = FALSE, verbose = TRUE, maxIter = params$bcvMaxIter))
          },
          message = function(m) {
            shinyjs::html(id = "bcvtable", html = tableHelper(m$message, nrow = 2),
                          add = FALSE)
            incProgress(1)
          })
          values$bcvRes <- res$counts
          values$bcvBest <- res$best
          values$bcvValid <- TRUE
        })
    }
  }
  
  observeEvent(
    input$bcvButton,
    {
      shinyjs::disable("bcvButton")
      shinyjs::disable("bcvAndBiclusterButton")
      guiBCV()
    }
  )
  
  observeEvent(input$bcvAndBiclusterButton, {
    shinyjs::disable("bcvButton")
    shinyjs::disable("bcvAndBiclusterButton")
    guiBCV()
    updateSliderInput(session, inputId = "k", value = as.integer(values$bcvBest))
    shinyjs::click("bicluster")
    updateNavbarPage(session, inputId = "navbarpage", 
                     selected = "Bicluster")
  })
  
  # A barplot of BCV results
  output$bcvPlot <- renderPlot({
    validate(need(length(values$bcvRes) > 0 && values$bcvBest > 0, ""))
    bcvPlotHelper()
  })
  bcvPlotHelper <- reactive({
    bcvRes <- values$bcvRes
    cols <- rep("#000000", length(bcvRes))
    cols[as.integer(values$bcvBest)] <- "#2ca25f"
    # assume that bcvRes is for k = 1:length(bcvRes
    barplot(bcvRes, xlab = "# Biclusters", ylab = "Times chosen as optimal", col = cols,
            axes = FALSE)
    axis(2, at = 0:4 * (ceiling(max(bcvRes) / 4)), 
         labels = as.integer(0:4 * (ceiling(max(bcvRes) / 4))))
  })
  
  observeEvent(values$bcvValid,
               {
                 if(values$bcvValid == FALSE) {
                   shinyjs::enable("bcvButton")
                   shinyjs::enable("bcvAndBiclusterButton")
                   shinyjs::show("bcvtable")
                   shinyjs::hide("bcvPlot")
                 } else {
                   shinyjs::hide("bcvtable")
                   shinyjs::show("bcvPlot")
                 }
               }
  )
  
  tableHelper <- function(str, nrow) {
    elements <- scan(text = str, quiet = TRUE)
    ncol <- ceiling(length(elements) / nrow)
    open <- "<table class=\"table table-hover\">\n<thead>\n"
    content <- do.call(paste0, lapply(seq_len(nrow), function(row) { # to produce each row
      paste0("<tr>\n", 
             do.call(paste0, lapply(seq_len(ncol), function(col) {
               # each th tag in the row
               paste0("<th>", elements[(row - 1) * ncol + col],  "</th>\n")
             })),
             # end of row or row header
             if(row == 1) "</tr>\n</thead>\n" else "</tr>\n"
      )
    }))
    close <- "</table>"
    paste0(open, content, close)
  }
  
  #### GO ENRICHMENT ####
  observeEvent(
    input$go,
    {
      shinyjs::disable("go")
      validate(need(inherits(values$bce, "BiclusterExperiment") &&
                      inherits(values$strategy, "BiclusterStrategy"),
                    ""))
      withProgress(
        message = "Searching for Gene Ontology enrichment...",
        value = 0, max = nclust(values$strategy),
        {
          values$goRes <- withCallingHandlers({
            testFE(bce = values$bce, strategy = values$strategy)
          },
          message = function(m) {
            incProgress(1)
            showNotification(conditionMessage(m), duration = 1)
          },
          warning = function(w) {
            showNotification(w$message, duration = NULL)
            },
          error = function(e) {
            showNotification(e$message, duration = NULL)
            return(NULL)}
          )
        })
      if(length(values$goRes) > 0) values$goValid = TRUE
    })
  
  observeEvent({
    values$goValid
    values$strategy
  }, { 
    if(values$goValid == FALSE && inherits(values$strategy, 
                                           "BiclusterStrategy")) {
      shinyjs::enable("go")
    } else {
      shinyjs::disable("go")
    }
  }
  )
  
  # a list of hypergeometric test result dataframes, one per bicluster
  goDF <- reactive({
    goRes <- values$goRes
    return(
      if(length(goRes) > 0) {
        lapply(goRes, function(listOfHyperGs) {
          p.value <- do.call(c, lapply(listOfHyperGs, GOstats::pvalues))
          adj.p.value <- p.adjust(p.value, method = "BH")
          odds.ratio <- do.call(c, lapply(listOfHyperGs, GOstats::oddsRatios))
          expected.genes <- do.call(c, lapply(listOfHyperGs, 
                                              GOstats::expectedCounts))
          matched.bicluster.ids <- do.call(c, lapply(listOfHyperGs, 
                                               Category::geneIdsByCategory))
          matched.dataset.ids <- do.call(c, lapply(listOfHyperGs,
                                                   Category::geneIdUniverse))
          matched.bicluster.genes <- unlist(lapply(matched.bicluster.ids,
                                                   length))
          matched.dataset.genes <- unlist(lapply(matched.dataset.ids, length))

          df <- data.frame(matched.bicluster.genes, matched.dataset.genes, 
                           expected.genes, odds.ratio, p.value, adj.p.value)
          df$matched.bicluster.ids <- matched.bicluster.ids
          df$matched.dataset.ids <- matched.dataset.ids
          df
        })
      } else list()
    )
  })
  
  goSummary <- reactive({
    goRes <- values$goRes
    if(length(goRes) > 1) {
      # summary data is same regardless of bicluster tested
      myHyperGResult <- goRes[[1]]
      # list of matched IDs for each GO term
      matchedDatasetIds <- Category::geneIdUniverse(myHyperGResult)
      # how many in the dataset were actually GO-annotated
      matchedDatasetSize <- length(unique(unlist(matched.dataset.ids)))
      
      # get the lowest adj. p-value for each bicluster
      significance <- sapply(goDF(), function(df) min(df$adj.p.value))
      if(any(significance < .Machine$double.eps)) {
        significance <- significance + .Machine$double.eps
      }
      significance <- sort(-log10(significance), decreasing = TRUE)
      
      return(list(universeSize = matchedDatasetSize,
                  significance = significance))
    }
  })
  
  # Plot of bicluster significance values
  output$goSigPlot <- renderPlot({
    goRes <- values$goRes
    validate(need(length(goRes) > 0 && !is.null(names(goRes)), ""))
    
    bp <- barplot(height = goSummary()$significance, xaxs = "i", yaxs = "i", xaxt = "n",
                  ylab = "-log10[p-value]")
    las <- if(any(nchar(names(goRes)) > 5)) 2 else 1 # labels rotated?
    axis(side = 1, at = bp, pos = 0, labels = names(goRes), las = las)
    abline(h = -log10(c(0.05)), col = "#2ca25f", lty = "dashed")
  })
  
  # Get input bicluster to display results for
  output$goBicluster <- renderUI({
    choices <- if(length(goDF()) > 0) names(goDF()) else list()
    val <- if(length(choices) > 0) choices[[1]] else NULL
    selectInput("goBicluster", label = "Select bicluster:",
                choices = choices, selected = val, multiple = FALSE)
  })
  
  # Display detailed GO results for one bicluster
  output$goTermTable <- DT::renderDT({
    validate(need(length(values$goRes) > 0 && !is.null(names(values$goRes)) &&
                    !is.null(input$goBicluster), "Please test GO enrichment"))
    df <- goDF()[[input$goBicluster]][, 1:6] # don't include gene lists
    return(
      DT::datatable(df, options = list(paging = FALSE), selection = 'single')
    )
  })
  
  observeEvent(
    input$goTabGenes,
    {selection <- input$goTermTable_rows_selected
    if(length(selection) > 0) {
      # assumes goDF() has not changed since goTermTable was rendered
      selection <- row.names(goDF()[[input$goBicluster]])[selection]
      
      updateSelectInput(session, inputId = "goTerm", selected = selection)
    }
    updateTabsetPanel(session, inputId = "goTab", selected = "Genes")
    })
  
  # Keep the choices in input$goTerm updated so that the "Inspect Genes"
  # button works.
  observeEvent({
    goDF()
    input$goBicluster
  }, {
    bicluster <- input$goBicluster
    if(!is.null(bicluster) && bicluster %in% names(goDF())) {
      choices <- row.names(goDF()[[bicluster]])
        updateSelectInput(session, inputId = "goTerm",
                          choices = choices,
                          selected = choices[1])
    }
  })
  
  output$goBiclusterGenes <- renderText({
    validate(need(length(goDF()) > 0, "Please test GO enrichment"))
    validate(need(!is.null(input$goTerm), "Please select a goTerm"))
    validate(need(input$goBicluster %in% names(goDF()), 
                  "Please select a valid bicluster"))
    
    genes <- goDF()[[input$goBicluster]][input$goTerm, "matched.bicluster.ids"]
    paste(unlist(genes), collapse = "\n")
  })
  
  output$goUniverseGenes <- renderText({
    validate(need(length(goDF()) > 0, "Please test GO enrichment"))
    validate(need(!is.null(input$goTerm), "Please select a goTerm"))
    validate(need(input$goBicluster %in% names(goDF()), 
                  "Please select a valid bicluster"))
    
    genes <- goDF()[[input$goBicluster]][input$goTerm, "matched.dataset.ids"]
    paste(unlist(genes), collapse = "\n")
  })
  
      # PANEL DESCRIPTIONS
      # output$explanation <- renderUI({
      #   res <- ""
      #   if (length(input$main_panel) > 0) {
      #     if ("Abundance" == input$main_panel) {
      #       res <- HTML(
      #         "Abundance shows a Z-score heatmap of the original matrix
      #         of metabolite peaks. Default ordering is distance-based by
      #         cluster membership. Sidebar options allow reordering
      #         of rows to highlight samples that are members of the
      #         currently selected cluster.
      #         Currently the annotations 'species' and 'tissue' are fictitious."
      #       )
      #     }
      # if ("Stability" == input$main_panel) {
      #   res <- HTML(
      #     "Stability index shows how stable each cluster
      #     is accross the selected range of <i>k</i>s.
      #     The stability index varies between 0 and 1, where
      #     1 means that the same cluster appears in every
      #     solution for different <i>k</i>."
      #   )
      # }
      # if ("Biomarkers" == input$main_panel) {
      #   res <- HTML(
      #     paste0(
      #       "This is a plot of loading, a statistic that
      #       quantifies the importance of each feature in distinguishing
      #       samples that are members of this bicluster."
      #     )
      #     )
      # }
      # if ("Bicluster members" == input$main_panel) {
      #   res <- HTML(
      #     "This is a score thresholding plot, showing which
      #     samples are above the threshold. The user should be
      #     told which thresholding method was used."
      #   )
      # }
      # return(res)
      # }
      # })
      
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