function(input, output, session) {
    params <- list(
        biclusterargs = if (dbg) {
            list(
                maxIter = 10L,
                withinVar = 10,
                row.release = 0.1,
                col.release = 0.1
            )
        } else
            list(),
        # currently otsu is hardcoded even in debug mode
        annotateBiclusters = if (dbg)
            3L
        else
            NA # how many biclusters to annotate
    )

    #### REACTIVE VALUES ####
    dep <- reactiveValues(
        pheatmap = requireNamespace("pheatmap", quietly = TRUE),
        BiocManager = requireNamespace("BiocManager", quietly = TRUE),
        gostats = requireNamespace("GOstats", quietly = TRUE)
    )

    values <- reactiveValues(
        # Initialize as the user's BiclusterExperiment.
        bce = if (inherits(userBce, "BiclusterExperiment")) {
            userBce
        } else
            NULL,
        strategy = if (inherits(userBce, "BiclusterExperiment")) {
            if (length(strategies(userBce)) > 0) {
                getStrat(userBce, 1)
            }
        } else
            NULL,
        zoom = NULL,
        bcvRes = numeric(0),
        bcvBest = 0,
        bcvValid = FALSE,
        goRes = numeric(0),
        # a named list containing results for each bicluster
        goLastParams = list()
    )

    kSliderProxy <- reactiveValues(value = NULL)
    observe({
        kSliderProxy$value <- input$k
    }, priority = 10)
    updateK <- function(value) {
        # updates the value of kSliderProxy immediately, without checking for
        # validity
        try({
            value <- as.integer(value)
            updateSliderInput(session, "k", value = value)
            kSliderProxy$value <- value
        })
    }

    output$status <- renderText({
        if (input$navbarpage == "Functional Annotation") {
            enrichedBc <- length(goSummary()$significance < 0.05)
            totalBc <- length(goSummary()$significance)
            if (inherits(goSummary()$universeSize, "integer")) {
                numMatched <- goSummary()$universeSize
            } else
                numMatched <- "--"
            desc <-
                paste(
                    enrichedBc,
                    "/",
                    totalBc,
                    "biclusters enriched (p < 0.05)",
                    "|",
                    numMatched,
                    "/",
                    nrow(values$bce),
                    "genes matched to GO terms"
                )
        } else if (input$navbarpage == "Data" ||
                   input$navbarpage == "Optimize") {
            if (inherits(values$bce, "BiclusterExperiment")) {
                samples <- ncol(values$bce)
                feat <- nrow(values$bce)
            } else {
                samples <- "--"
                feat <- "--"
            }
            desc <- paste(samples, "samples X", feat, "feature")
        } else if (input$navbarpage == "Bicluster") {
            if (inherits(values$strategy, "BiclusterStrategy")) {
                k <- nclust(values$strategy)
            } else {
                k <- "--"
            }
            desc <-
                paste(k, "Biclusters in the selected BiclusterStrategy")
        } else {
            desc <- ""
        }
        return(desc)
    })

    #### DATA I/O ####
    observeEvent({
        input$input_df
        input$sepchar
        input$decchar
        input$quotechar
        input$skiplines
        input$row_names
        input$header
    },
    {
        withProgress({
            incProgress(1 / 8)
            if (!is.null(input$input_df)) {
                rows <- if (input$row_names)
                    1
                else
                    NULL
                rawmat <-
                    read.table(
                        input$input_df$datapath,
                        header = input$header,
                        row.names = rows,
                        fill = TRUE,
                        comment.char = "",
                        sep = input$sepchar,
                        quote = input$quotechar,
                        dec = input$decchar,
                        skip = input$skiplines
                    )
                incProgress(4 / 8)
                rawmat <- as.matrix(rawmat)
                bce <- BiclusterExperiment(rawmat)
            } else {
                # no file -> load user BiclusterExperiment
                bce <- userBce
            }
            incProgress(1 / 8)
        }, message = "Parsing data...", value = 0)
        values$bce <- bce
        values$bcvValid <- FALSE
    })

    observeEvent(input$postUploadUpdate, {
        validate(need(inherits(
            values$bce, "BiclusterExperiment"
        ), ""))
        withProgress(message = "Updating", value = 1 / 8, {
            bce <- values$bce
            if (input$transpose) {
                # careful; is all data preserved here?
                bce <- BiclusterExperiment(t(as.matrix(userBce)))
            }
            setProgress(3 / 8)
            if (nchar(input$customRowNames) > 0) {
                ns <- unlist(strsplit(
                    x = input$customRowNames,
                    split = "\\s+"
                ))
                if (length(ns) == dim(bce)[1]) {
                    featureNames(bce) <- ns
                }
            }
            setProgress(5 / 8)
            if (nchar(input$customColNames) > 0) {
                ns <- unlist(strsplit(
                    x = input$customColNames,
                    split = "\\s+"
                ))
                if (length(ns) == dim(bce)[2]) {
                    sampleNames(bce) <- ns
                }
            }
            setProgress(7 / 8)
            values$bce <- bce
            values$bcvValid <- FALSE
        })
    })

    reactiveSeparator <- observeEvent(input$input_df, {
        # If the file input changes, try automatically changing separator to ","
        if (substr(
            input$input_df$name,
            nchar(input$input_df$name) - 2,
            nchar(input$input_df$name)
        ) == "csv") {
            updateTextInput(session, "sepchar", value = ",")
        }
    })

    # Update this central matrix from its source
    controls <- reactive({
        if (is.null(input$input_df)) {
            # If the user has not chosen a file yet, look for a
            # BiclusterExperiment in the execution environment of biclusterGUI
            shinyjs::disable("row_names")
            shinyjs::disable("header")
            shinyjs::disable("skiplines")
            shinyjs::disable("sepchar")
            shinyjs::disable("quotechar")
            shinyjs::disable("decchar")
        } else {
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
        }
    })

    # Render an interactive data.table for the user's benefit
    output$dt <- suppressWarnings(
        DT::renderDT({
            validate(need(
                inherits(values$bce, "BiclusterExperiment"),
                "You may import your dataset"
            ))
            return(as.matrix(values$bce))
        }, options = list(paging = FALSE, info = FALSE), fillContainer = TRUE,
        autoHideNavigation = TRUE, server = TRUE, selection = "none")
    )

    # helper functions
    reactiveHeatmapHeight500 <- reactive({
        max(500, length(input$annots) * 33 +
                sum(unlist(
                    lapply(input$annots, function(annot) {
                        length(unique(cbind(
                            Biobase::pData(Biobase::phenoData(values$bce)),
                            pred(getStrat(values$bce, input$strategy))
                        )[, annot]))
                    })
                )) * 22 - 106)
    })

    output$annotPicker <- renderUI({
        classes <- if (inherits(values$bce, "BiclusterExperiment")) {
            c(colnames(Biobase::phenoData(values$bce)))
        } else
            NULL
        selectInput(
            "annots",
            label = "Annotations",
            choices = classes,
            selected = NULL,
            multiple = TRUE
        )
    })

    # plot abundance heatmap (original data)
    output$uiabundance <- renderUI({
        validate(need(
            inherits(values$bce, "BiclusterExperiment"),
            "You may import your dataset."
        ))
        # withProgress(message = "Plotting...", value = 3/8, {
        height = reactiveHeatmapHeight500()
        return(plotOutput("abundance", height = height))
        # })
    })
    output$abundance <- renderPlot({
        gt <- reactive_abundance()
        print(gt) # printing ensures the returned gTable is drawn to Shiny
    }, height = function() {
        reactiveHeatmapHeight500()
    })
    reactive_abundance <- reactive({
        bce <- values$bce

        if (input$logBase == "e") {
            logBase <- exp(1)
        } else {
            logBase <- as.numeric(input$logBase)
        }
        phenoLabels <-
            intersect(input$annots, colnames(Biobase::phenoData(bce)))
        # biclustLabels <- if(length(names(bce)) > 0) {
        #   intersect(input$annots, bcNames(getStrat(bce, input$strategy)))
        plot(
            bce,
            logBase = logBase,
            phenoLabels = phenoLabels,
            ordering = if (input$heatmapReorder)
                "distance"
            else
                "input",
            strategy = input$strategy,
            rowNames = input$featNames,
            colNames = input$sampNames
        )
    })

    # Plot of samples along first two PCs
    output$pca <- renderPlot({
        validate(need(
            inherits(values$bce, "BiclusterExperiment"),
            "You may import your dataset."
        ))
        withProgress(message = "Plotting...", value = 3 / 8, {
            pca(values$bce)
        })
    })

    # Euclidean Distance between samples (applied to raw data...might be good to
    # scale first?)
    output$sampleDistance <- renderPlot({
        validate(need(
            inherits(values$bce, "BiclusterExperiment"),
            "You may import your dataset."
        ))
        gt <- sampleDistance()
        print(gt)
    }, height = function() {
        reactiveHeatmapHeight500()
    })
    sampleDistance <- reactive({
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
            phenoLabels <-
                intersect(input$annots, colnames(Biobase::phenoData(bce)))
            distType <-
                if (input$sampDistType == "Euclidean distance") {
                    "euclidean"
                } else {
                    "pearson"
                }
            plotDist(
                bce,
                type = "sample",
                distType = distType,
                phenoLabels = phenoLabels,
                biclustLabels = NULL,
                ordering = if (input$sampDistReorder)
                    "distance"
                else
                    "input",
                strategy = input$strategy,
                rowColNames = input$sampDistRownames
            )
        })
    })

    #### BICLUSTER TAB ####
    # FIXME BCV and Bicluster still cannot modify the sliderInput.
    observeEvent({
        values$bce
    },
    {
        bce <- values$bce
        val <- kSliderProxy$value
        if (is.null(bce)) {
            maxk <- 1
            val <- 1
        } else {
            m <- as.matrix(bce)
            maxk <- min(nrow(m), ncol(m))
            val <- min(2, maxk) # if possible, suggest 2 to the user
        }
        if (inherits(values$strategy, "BiclusterStrategy")) {
            updateSelectInput(session, "algo",
                              selected = method(values$strategy))
            # if there's already a biclusterstrategy..
            val <- nclust(values$strategy)
        }
        updateSliderInput(session, "k", max = maxk, step = 1)
        updateK(value = val)
    })

    observeEvent({
        values$bce
        kSliderProxy$value
        input$algo
    },
    {
        if (biclusteringAllowed()) {
            shinyjs::enable("bicluster")
        } else {
            shinyjs::disable("bicluster")
        }
    })
    biclusteringAllowed <- function() {
        # core data must be available
        if (inherits(values$bce, "BiclusterExperiment") &&
            inherits(kSliderProxy$value, "integer")) {
            # The parameters must not match an existing BiclusterStrategy
            if (all(names(values$bce) !=
                    name(
                        list(
                            bca = capitalize(input$algo),
                            threshAlgo = capitalize("otsu"),
                            k = kSliderProxy$value
                        )
                    ))) {
                return(TRUE)
            }
        }
        return(FALSE)
    }

    observeEvent(input$bicluster, {
        validate(need(
            inherits(values$bce, "BiclusterExperiment"),
            "You may import your dataset."
        ))
        validate(need(inherits(kSliderProxy$value, "integer"), ""))
        # temporarily override the button's own renderUI
        shinyjs::disable("bicluster")
        try(removeNotification("overlapNotif"))

        bce <- values$bce
        # Look for an existing BiclusterStrategy with the same parameters
        stratName <- name(
            list(
                bca = capitalize(input$algo),
                sta = capitalize("otsu"),
                lta = capitalize("otsu"),
                k = kSliderProxy$value
            )
        )
        matchStrats <- which(names(bce) == stratName)
        if (length(matchStrats) == 0) {
            # if the requested BiclusterStrategy doesn't exist yet, run
            # addStrat()
            withProgress({
                # Hack to set withinVar dynamically
                if ("withinVar" %in% names(params$biclusterargs)) {
                    params$biclusterargs$withinVar <-
                        params$biclusterargs$withinVar * nrow(bce)
                }
                # append any optional debug-mode arguments
                withCallingHandlers({
                    set.seed(12345)
                    bce <- do.call(
                        addStrat,
                        c(
                            bce = bce,
                            k = kSliderProxy$value,
                            method = input$algo,
                            verbose = FALSE,
                            params$biclusterargs
                        )
                    )},
                    warning = function(w) {
                        # In case less than the requested number of biclusters
                        # was found
                        showNotification(w$message, duration = NULL)
                    }
                )
                newStrat <- names(bce)[length(names(bce))]
                algo <- strsplit(newStrat, split = " | ")[[1]][1]

                if (algo != capitalize(input$algo)) {
                    updateSelectInput(session,
                                      inputId = "algo",
                                      selected = algo)
                    showNotification(
                        paste(
                            input$algo,
                            "failed on your dataset, so the",
                            algo,
                            "algorithm was used instead."
                        ),
                        duration = NULL
                    )
                }
                # Whatever the new BiclusterStrategy is, set it as the active
                # strategy
                values$strategy <- getStrat(bce, newStrat)
                values$bce <- bce
            }, message = "Biclustering...", value = 0.2)
        } else {
            values$strategy <- getStrat(bce, matchStrats[1])
        }
    })

    # When the parameters are changed, see if biclustering has already been run
    # yet
    observeEvent({
        values$bce
        input$algo
        kSliderProxy$value
    }, {
        if (inherits(values$bce, "BiclusterExperiment") &&
            inherits(kSliderProxy$value, "integer")) {
            bce <- values$bce
            stratName <- name(list(
                bca = capitalize(input$algo),
                threshAlgo = capitalize("otsu"),
                k = kSliderProxy$value
            ))
            matchStrats <- which(names(bce) == stratName)
            if (length(matchStrats) > 0) {
                # If the BiclusterStrategy has been run already, activate it
                values$strategy <- strategies(bce)[[stratName]]
            } else {
                # If not, then erase the active strategy
                values$strategy <- NULL
            }
        }
    })

    #### Summary ####
    # Reset zoom when the user's input data changes or upon double-click
    observeEvent({
        values$bce
        input$bcHighlights_dblclick
    }, {
        validate(need(inherits(
            values$bce, "BiclusterExperiment"
        ), ""))
        mat <- as.matrix(values$bce)

        if (inherits(mat, "matrix") && mode(mat) == "numeric") {
            arr <- imageArr()
            values$zoom <- c(1, dim(arr)[2], 1, dim(arr)[1])
        }
    })

    # Re-render the bicluster plot when raw data OR BiclusterStrategy OR
    # selected biclusters changes
    imageArr <- reactive({
        # render the whole heatmap
        validate(need(inherits(
            values$bce, "BiclusterExperiment"
        ), ""))
        mat <- as.matrix(values$bce)
        validate(need(inherits(mat, "matrix") &&
                          mode(mat) == "numeric", ""))

        m <- mat / max(mat)
        width <- ncol(m)
        height <- nrow(m)
        # Convert the vector to an array with 3 planes for R, G, B
        arr <- array(c(m, m, m), dim = c(height, width, 3))

        # Add bicluster outlines
        bcs <- values$strategy
        if (inherits(bcs, "BiclusterStrategy")) {
            cols <- hcl(
                h = seq(0, (nclust(bcs) - 1) / (nclust(bcs)),
                        length = nclust(bcs)) * 360,
                c = 100,
                l = 65,
                fixup = TRUE
            )
            cols <-
                col2rgb(cols) # RGB as rows; each column a different color
            # FIXME: Allow to select a subset of biclusters. Allow bicluster ID
            # when hovering mouse over.
            biclusterLists <- biclusterMatrix2List(
                rowxBicluster = clusteredFeatures(bcs),
                biclusterxCol = clusteredSamples(bcs))
            displayOrder <-
                order(
                    sizes(
                        biclusterRows = biclusterLists[[1]],
                        biclusterCols = biclusterLists[[2]]
                    ),
                    decreasing = TRUE
                )
            lapply(displayOrder, function(bicluster) {
                yrange <- which(clusteredSamples(bcs)[bicluster,])
                xrange <- which(clusteredFeatures(bcs)[, bicluster])
                arr[xrange, yrange, 1] <<-
                    cols["red", bicluster] / 255
                arr[xrange, yrange, 2] <<-
                    cols["green", bicluster] / 255
                arr[xrange, yrange, 3] <<-
                    cols["blue", bicluster] / 255
            })
        }
        arr
    })

    output$bcHighlights <- renderImage({
        arr <- imageArr()

        # These may be specified in biclusterUI.R
        width <- session$clientData$output_bcHighlights_width
        height <- session$clientData$output_bcHighlights_height

        validate(
            need(
                inherits(arr, "array") &&
                    mode(arr) == "numeric" &&
                    !is.null(dim(arr)),
                "Data for heatmap has not been parsed yet or contains non-numbers"
            )
        )
        temp <- tempfile(fileext = ".png")

        # subset the PNG
        if(is.null(values$zoom)) values$zoom <- c(1, dim(arr)[2],
                                                  1, dim(arr)[1])

        rangeY <-
            c(values$zoom[3], values$zoom[4])
        rangeX <-
            c(values$zoom[1], values$zoom[2])

        arr <- arr[rangeY[1]:rangeY[2], rangeX[1]:rangeX[2], , drop = FALSE]

        # resize the subsetted array to the desired display resolution
        arr <- EBImage::resize(x = arr, w = height, h = width,
                               antialias = FALSE, filter = "none")
        try(png::writePNG(arr, target = temp))

        height <- input$dimension[2]
        # Return a list containing information about the zoomed image
        list(
            src = temp,
            contentType = "image/png",
            title = "Click and drag to zoom in; double-click to reset"
        )
    })

    # If the user uses the brush, zoom in. Clearing the brush simply allows the
    # GUI to await another brush input.
    observeEvent(input$bcHighlights_brush, {
        brush <- input$bcHighlights_brush
        zoom <- values$zoom

        if (!is.null(brush)) {
            # contains fractions of the whole image
            values$zoom <-
                c(
                    zoom[1] + floor((brush$xmin / brush$domain$right) *
                                        (zoom[2] - zoom[1] + 1)),
                    zoom[1] - 1 + ceiling((brush$xmax / brush$domain$right) *
                                              (zoom[2] - zoom[1] + 1)),
                    zoom[3] + floor((brush$ymin / brush$domain$bottom) *
                                        (zoom[4] - zoom[3] + 1)),
                    zoom[3] - 1 + ceiling((brush$ymax / brush$domain$bottom) *
                                              (zoom[4] - zoom[3] + 1))
                )
            shinyjs::runjs("document.getElementById('bcHighlights_brush').remove()")
        }
    })

    # Warn if biclusters overlap
    overlapWarn <- observeEvent(values$strategy, {
        validate(need(inherits(
            values$strategy, "BiclusterStrategy"
        ), ""))
        bcs <- values$strategy
        # Create lists of rows and columns contained in biclusters
        if (nclust(bcs) > 1) {

            rc <- biclusterMatrix2List(clusteredFeatures(bcs),
                                       clusteredSamples(bcs))
            biclusterRows <- rc[[1]]
            biclusterCols <- rc[[2]]
            # Check for overlaps with any other biclusters. Since a bicluster
            # always overlaps with itself completely, any overlap causes sum(l)
            # > 1
            overlaps <- overlap(biclusterRows, biclusterCols)
            nonOverlap <-
                vapply(
                    overlaps,
                    FUN.VALUE = numeric(1),
                    FUN = function(l)
                        sum(l == 1)
                )
            if (any(!nonOverlap)) {
                showNotification(
                    id = "overlapNotif",
                    paste(
                        do.call(paste, c(
                            as.list(bcNames(bcs)[!nonOverlap]), sep = ", "
                        )),
                        "overlap. Please interpret the heatmap annotations with care."
                    ),
                    duration = NULL
                )
            }
        }
    })

    # Download a BiclusterExperiment with all BiclusterStrategy objects
    output$downloadBceAllButton <- renderUI({
        if (inherits(values$bce, "BiclusterExperiment")) {
            if (length(strategies(values$bce)) > 0) {
                return(
                    downloadButton(
                        'downloadBceAll',
                        'Download all bicluster results',
                        class = "dlButton"
                    )
                )
            }
        }
        return(
            actionButton(
                'downloadBceAll',
                'Download all bicluster results',
                class = "dlButton",
                disabled = TRUE
            )
        )
    })
    output$downloadBceAll <- downloadHandler(
        filename = function() {
            paste0(substr(Sys.time(), start = 1, stop = 10),
                   "_BiclusterExperiment.Rdata")
        },
        content = function(file) {
            bce <- values$bce
            save(bce, file = file)
        }
    )

    # Download a BiclusterExperiment with the current BiclusterStrategy
    output$downloadBceButton <- renderUI({
        if (inherits(values$strategy, "BiclusterStrategy")) {
            downloadButton('downloadBceCurrent',
                           'Download this bicluster result',
                           class = "dlButton")
        } else {
            actionButton(
                'downloadBceCurrent',
                'Download this bicluster result',
                class = "dlButton",
                disabled = TRUE
            )
        }
    })
    output$downloadBceCurrent <- downloadHandler(
        filename = function() {
            paste0(substr(Sys.time(), start = 1, stop = 10),
                   "_BiclusterExperiment.Rdata")
        },
        content = function(file) {
            bce <- wipeExcept(values$bce, values$strategy)
            save(bce, file = file)
        }
    )

    #### Samples ####

    #would be nice to Invalidate these panels when user changes kSlider, but
    #hasn't clicked "Run" yet
    output$sampleBicluster <- renderUI({
        choices <- if (inherits(values$strategy, "BiclusterStrategy")) {
            bcNames(values$strategy)
        } else {
            list()
        }
        selectInput("sampleBicluster", "Select bicluster:",
                    choices = choices)
    })

    # Heatmap of loadings for all samples
    output$sampleHeatmap <- renderPlot({
        gt <- reactivesampleHeatmap()
        print(gt)
    })
    reactivesampleHeatmap <- reactive({
        validate(need(
            inherits(values$strategy, "BiclusterStrategy"),
            "Please bicluster first."
        ))
        withProgress(message = "Plotting...",
                     value = 0, {
                         set.seed(1234567)
                         factorHeatmap(
                             bce = values$bce,
                             bcs = values$strategy,
                             type = "sample",
                             ordering = if (input$sampleReorder) {
                                 "cluster"
                             } else {
                                 "input"
                             },
                             colNames = input$biclusterFeatNames
                         )
                     })
    })

    # Plot of sample loadings for one bicluster
    output$samplePlot <- renderPlot({
        validate(need(
            inherits(values$strategy, "BiclusterStrategy") &&
                nchar(input$sampleBicluster) > 0,
            ""
        ))
        reactiveMarkers()
    })
    reactiveMarkers <- reactive({
        withProgress(message = "Plotting...", value = 0, {
            set.seed(1234567)
            plotThreshold(
                bce = values$bce,
                bcs = values$strategy,
                type = "sample",
                bicluster = input$sampleBicluster,
                ordering = if (input$sampleReorder) {
                    "cluster"
                } else {
                    "input"
                },
                xlabs = input$biclusterFeatNames
            )
        })
    })

    # Get sample list for selected bicluster
    output$biclusterSampleList <- renderText({
        bce <- values$bce
        bcs <-  values$strategy
        bicluster <- input$sampleBicluster
        validate(need(
            inherits(bce, "BiclusterExperiment") &&
                !is.null(bicluster) &&
                inherits(bcs, "BiclusterStrategy"),
            ""
        ))
        validate(need(nchar(bicluster) > 0, ""))

        sampleI <- which(clusteredSamples(bcs)[bicluster,])
        samples <- sampleNames(bce)[sampleI]
        paste(unlist(samples), collapse = "\n")
    })

    output$biclusterSampleListLabel <- renderUI({
        bicluster <- input$sampleBicluster
        desc <- if (is.null(bicluster)) {
            "Biclustered samples:"
        } else if (nchar(bicluster) == 0) {
            "Biclustered samples:"
        } else {
            paste(bicluster, "samples:")
        }
        tags$h4(desc)
    })

    #### Features ####
    output$featureBicluster <- renderUI({
        choices <- if (inherits(values$strategy, "BiclusterStrategy")) {
            bcNames(values$strategy)
        } else {
            list()
        }
        selectInput("featureBicluster", "Select bicluster:",
                    choices = choices)
    })

    reactiveHeatmapHeight300 <- reactive({
        bce <- values$bce
        max(300, length(input$annots) * 33 +
                sum(unlist(
                    lapply(input$annots, function(annot) {
                        length(unique(cbind(
                            Biobase::pData(Biobase::phenoData(bce)),
                            pred(getStrat(bce, input$strategy))
                        )[, annot]))
                    })
                )) * 22 - 106)
    })

    # heatmap of scores for all features
    # using renderUI prevents overlapping plots
    output$featureHeatmap <- renderUI({
        height = reactiveHeatmapHeight300()
        plotOutput("featureHeatmap", height = height)
    })
    output$featureHeatmap <- renderPlot({
        gt <- reactivefeatureHeatmap()
        print(gt)
    },
    height = function() {
        reactiveHeatmapHeight300()
    })
    reactivefeatureHeatmap <- reactive({
        validate(need(
            inherits(values$strategy, "BiclusterStrategy"),
            "Please run biclustering first."
        ))
        bce <- values$bce
        withProgress(message = "Plotting...",
                     value = 0, {
                         set.seed(1234567)
                         # Create the annotation selector
                         phenoLabels <- intersect(input$annots,
                                                  colnames(Biobase::phenoData(
                                                      bce)))
                         biclustLabels <- intersect(input$annots,
                                                    bcNames(values$strategy))
                         factorHeatmap(
                             bce = bce,
                             bcs = values$strategy,
                             type = "feature",
                             phenoLabels = phenoLabels,
                             biclustLabels = biclustLabels,
                             ordering = if (input$featureReorder) {
                                 "cluster"
                             } else {
                                 "input"
                             },
                             colNames = input$biclusterSampNames
                         )
                     })
    })

    # Plot of feature scores for one bicluster
    output$featurePlot <- renderPlot({
        featurePlotHelper()
    })
    featurePlotHelper <- reactive({
        validate(need(
            inherits(values$strategy, "BiclusterStrategy") &&
                !is.null(input$featureBicluster),
            ""
        ))
        withProgress(message = "Plotting...", value = 0, {
            set.seed(1234567)
            plotThreshold(
                bce = values$bce,
                bcs = values$strategy,
                type = "feature",
                bicluster = input$featureBicluster,
                ordering = if (input$featureReorder) {
                    "cluster"
                } else {
                    "input"
                },
                xlabs = input$biclusterSampNames
            )
        })
    })

    # Get marker list for selected bicluster
    output$biclusterFeatureList <- renderText({
        bce <- values$bce
        bcs <-  values$strategy
        bicluster <- input$featureBicluster
        validate(need(
            inherits(bce, "BiclusterExperiment") &&
                !is.null(bicluster) &&
                inherits(bcs, "BiclusterStrategy"),
            ""
        ))
        validate(need(nchar(bicluster) > 0, ""))
        geneI <- which(clusteredFeatures(bcs)[, bicluster])
        genes <- featureNames(bce)[geneI]
        paste(unlist(genes), collapse = "\n")
    })

    output$biclusterFeatureListLabel <- renderUI({
        bicluster <- input$featureBicluster
        desc <- if (is.null(bicluster)) {
            "Biclustered features:"
        } else if (nchar(bicluster) == 0) {
            "Biclustered features:"
        } else {
            paste(bicluster, "markers:")
        }
        tags$h4(desc)
    })

    #### BI-CROSS-VALIDATION ####
    output$bcvButton <- renderUI({
        if (!inherits(values$bce, "BiclusterExperiment") ||
            values$bcvValid) {
            # disabled button if invalid data, or BCV already run
            return(actionButton(
                inputId = "bcvButton",
                label = "Perform BCV",
                disabled = TRUE
            ))
        } else if (any(is.na(as.matrix(values$bce)))) {
            # button redirecting to confirmation dialog
            return(actionButton(inputId = "bcvCheckButton",
                                label = "Perform BCV"))
        } else {
            # BCV immediately
            return(actionButton(inputId = "bcvButton", label = "Perform BCV"))
        }
    })
    output$bcvAndBiclusterButton <- renderUI({
        if (!inherits(values$bce, "BiclusterExperiment") ||
            values$bcvValid) {
            return(
                actionButton(
                    inputId = "bcvAndBicluster",
                    label = "Perform BCV and bicluster",
                    disabled = TRUE
                )
            )
        } else if (any(is.na(as.matrix(values$bce)))) {
            return(
                actionButton(inputId = "bcvCheckAndBicluster",
                             label = "Perform BCV and bicluster")
            )
        } else {
            return(
                actionButton(inputId = "bcvAndBicluster",
                             label = "Perform BCV and bicluster")
            )
        }
    })

    observeEvent(input$bcvCheckButton, {
        try(removeNotification("bcvNaNotif"))
        showNotification(
            id = "bcvNaNotif",
            ui = paste(
                "Since some matrix elements are NA, bi-cross-validation",
                "will run\n much more slowly."
            ),
            action = actionButton(inputId = "bcvButton", label = "Continue"),
            duration = NULL
        )
    })
    observeEvent(input$bcvButton,
                 {
                     runBcv()
                     updateK(val = values$bcvBest)
                 })
    runBcv <- function() {
        shinyjs::disable("bcvButton")
        shinyjs::disable("bcvAndBiclusterButton")
        validate(need(inherits(
            values$bce, "BiclusterExperiment"
        ), ""))
        withProgress(
            message = "Performing bi-cross-validation...",
            value = 0,
            max = params$bcvMaxIter,
            {
                res <- withCallingHandlers({
                    shinyjs::html("bcvtext", "")
                    suppressWarnings(do.call(auto_bcv,
                                             c(
                                                 list(
                                                     Y = as.matrix(values$bce),
                                                     bestOnly = FALSE,
                                                     verbose = TRUE,
                                                     interactive = FALSE
                                                 )
                                             )))
                },
                message = function(m) {
                    shinyjs::html(
                        id = "bcvtable",
                        html = tableHelper(m$message, nrow = 2),
                        add = FALSE
                    )
                    incProgress(1)
                })
                values$bcvRes <- res$counts
                values$bcvBest <- res$best
                values$bcvValid <- TRUE
            }
        )
    }

    observeEvent(input$bcvCheckAndBicluster,
                 {
                     try(removeNotification("bcvNaNotifThenBicluster"))
                     showNotification(
                         id = "bcvNaNotifThenBicluster",
                         ui = paste(
                             "Since some matrix elements are NA, ",
                             "bi-cross-validation will run\n much more slowly."
                         ),
                         action = actionButton(inputId = "bcvAndBicluster",
                                               label = "Continue"),
                         duration = NULL
                     )
                 })
    observeEvent(input$bcvAndBicluster, {
        runBcv()

        updateK(value = values$bcvBest)
        if (biclusteringAllowed()) {
            shinyjs::enable("bicluster")
            shinyjs::click("bicluster")
        }
        updateTabsetPanel(session, "navbarpage", selected = "Bicluster")
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
        barplot(
            bcvRes,
            xlab = "# Biclusters",
            ylab = "Times chosen as optimal",
            col = cols,
            axes = FALSE
        )
        axis(2,
             at = 0:4 * (ceiling(max(bcvRes) / 4)),
             labels = as.integer(0:4 * (ceiling(
                 max(bcvRes) / 4
             ))))
    })

    observeEvent(values$bcvValid,
                 {
                     if (values$bcvValid == FALSE) {
                         shinyjs::enable("bcvButton")
                         shinyjs::enable("bcvAndBiclusterButton")
                         shinyjs::show("bcvtable")
                         shinyjs::hide("bcvPlot")
                     } else {
                         shinyjs::hide("bcvtable")
                         shinyjs::show("bcvPlot")
                     }
                 })

    tableHelper <- function(str, nrow) {
        elements <- scan(text = str, quiet = TRUE)
        ncol <- ceiling(length(elements) / nrow)
        open <- "<table class=\"table table-hover\">\n<thead>\n"
        content <-
            do.call(paste0, lapply(seq_len(nrow), function(row) {
                # to produce each row
                paste0(
                    "<tr>\n",
                    do.call(paste0, lapply(seq_len(ncol), function(col) {
                        # each th tag in the row
                        paste0("<th>", elements[(row - 1) * ncol + col],
                               "</th>\n")
                    })),
                    # end of row or row header
                    if (row == 1)
                        "</tr>\n</thead>\n"
                    else
                        "</tr>\n"
                )
            }))
        close <- "</table>"
        paste0(open, content, close)
    }

    #### GO ENRICHMENT ####
    # FIXME Don't crash if gene IDs aren't ENSEMBL
    output$go <- renderUI({
        if (!dep$gostats) {
            return(actionButton("gostats", "Install dependency \"GOstats\""))
        } else {
            # all input must be present, and params must have changed since last
            # run
            active <- goAllowed() &&
                !identical(list(values$strategy, input$orgDb, input$gos),
                           values$goLastParams)
            if (active) {
                return(actionButton("go", "Test for GO enrichment"))
            } else {
                return(actionButton("go", "Test for GO enrichment",
                                    disabled = TRUE))
            }
        }
    })
    goAllowed <- reactive({
        if (inherits(values$strategy, "BiclusterStrategy")) {
            allowed <- TRUE
            # the species dropdown must have loaded
            if (is.null(input$orgDb)) {
                allowed <- FALSE
            } else {
                if (nchar(input$orgDb) == 0) {
                    allowed <- FALSE
                }
            }
        } else {
            # there must be biclustering data available to get gene lists from
            allowed <- FALSE
        }
        return(allowed)
    })
    # Call this anytime BiocManager is absent
    requestBiocManager <- function() {
        if (!dep$BiocManager) {
            showNotification(
                ui = paste(
                    "BiocManager must be installed. Clicking",
                    "below will run",
                    "install.packages(\"BiocManager\", quiet = TRUE)"
                ),
                action = actionButton("biocmanager", "Install"),
                id = "biocManagerNotif",
                duration = NULL,
                closeButton = FALSE
            )
            # I hate to make the notification non-closable, but I have no way of
            # detecting the close action. Then the user has no way to install
            # BiocManager after closing.
        }
    }
    # Install GOstats
    observeEvent(input$gostats,
                 {
                     if (!dep$BiocManager) {
                         requestBiocManager()
                     } else {
                         shinyjs::disable("gostats") # user may click only once
                         withProgress(
                             message = "Installing GOstats...",
                             value = 3 / 8, {
                                 suppressWarnings(BiocManager::install(
                                     "GOstats", update = FALSE, ask = FALSE))
                             })
                         withProgress(
                             message = "Verifying installation...",
                             value = 7 / 8, {
                                 if (requireNamespace("GOstats",
                                                      quietly = TRUE))
                                 {
                                     # trigger update of the "Run" button
                                     dep$gostats <- TRUE
                                     showNotification( ui = "GOstats installed",
                                                       duration = NULL)
                                 } else {
                                     showNotification(
                                         ui = paste(
                                             "Unable to install GOstats. ",
                                             "Please save your work exit the",
                                             " GUI, and install GOstats",
                                             " manually before attempting to",
                                             " test for functional enrichment."
                                         ), duration = NULL)
                                 }
                             })
                     }
                     # cleanup? will this button ever display again?
                     shinyjs::enable("gostats")
                 })
    observeEvent(
        input$biocmanager,
        {
            removeNotification("biocManagerNotif")
            withProgress(
                message = "Running install.packages(\"BiocManager\", quiet = TRUE)",
                value = 3 / 8,
                {
                    suppressWarnings(install.packages("BiocManager",
                                                      quiet = TRUE))
                })
            if (requireNamespace("BiocManager", quietly = TRUE)) {
                dep$BiocManager <- TRUE
            }
        })
    observeEvent(
        input$go,
        {
            validate(need(
                inherits(values$bce, "BiclusterExperiment") &&
                    inherits(values$strategy, "BiclusterStrategy") &&
                    length(input$gos) > 0,
                ""
            ))

            # temporarily override the button's own renderUI
            shinyjs::disable("go")

            if (!requireNamespace(input$orgDb, quietly = TRUE)) {
                try(removeNotification("orgDbInstallNotif"))

                # If the database needs to be installed, help the user to do so.
                showNotification(
                    paste(input$orgDb, "must be installed from Bioconductor"),
                    action = actionButton("orgDbInstall", "Install"),
                    id = "orgDbInstallNotif",
                    duration = NULL
                )

                # the user can click this button repeatedly, but only one notif
                # is shown
                shinyjs::enable("go")
                # User will have to click the button again AFTER installation

            } else {
                withProgress(
                    message = "Searching for Gene Ontology enrichment...",
                    value = 0,
                    max = nclust(values$strategy) * length(input$gos),
                    {
                        values$goRes <- withCallingHandlers({
                            # withCallingHandlers handles warnings; errors are
                            # trapped by try()
                            try(testFE(
                                bce = values$bce,
                                strategy = values$strategy,
                                go = input$gos,
                                orgDb = input$orgDb
                            ))
                        },
                        message = function(m) {
                            incProgress(1)
                            showNotification(conditionMessage(m), duration = 5)
                        },
                        warning = function(w) {
                        })
                        if (inherits(values$goRes, "try-error")) {
                            showNotification(values$goRes, duration = NULL)
                            values$goRes <-
                                NULL # return to previous app state
                        } else {
                            # store the running parameters
                            values$goLastParams <-
                                list(values$strategy, input$orgDb, input$gos)
                        }
                    }
                )
            }
        })
    output$species <- renderUI({
        if (!dep$BiocManager) {
            # pop up requesting BiocManager. Selector will be empty until
            # installed.
            requestBiocManager()
            return(selectInput(
                "orgDb",
                "Species:",
                choices = NULL,
                selected = NULL
            ))
        } else {
            # search for Org.*.db packages on Bioconductor
            pkgs <-
                available.packages(
                    repo = BiocManager::repositories()["BioCann"],
                    type = "source")
            orgdbI <- grep(pattern = "^org.*", row.names(pkgs))
            return(selectInput("orgDb", "Org. database:", choices =
                                   row.names(pkgs[orgdbI,])))
        }
    })
    # Install the selected org.*.Db
    observeEvent(
        input$orgDbInstall,
        {
            if (dep$BiocManager) {
                # this check is redundant
                try(removeNotification("orgDbInstallNotif"))
                withProgress(message = paste("Installing", input$orgDb),
                             value = 3 / 8,
                             {
                                 BiocManager::install(
                                     input$orgDb,
                                     update = FALSE,
                                     ask = FALSE,
                                     suppressAutoupdate = TRUE,
                                     type = "source"
                                 )
                             })
                if (requireNamespace(input$orgDb, quietly = TRUE)) {
                    # why doesn't this work? shinyjs::click("go")
                    showNotification(paste(input$orgDb, "installed!"),
                                     duration = NULL)
                }
            } else {
                requestBiocManager()
            }
        })

    # a list of hypergeometric test result dataframes, one per bicluster
    goDF <- reactive({
        goRes <- values$goRes
        return(if (length(goRes) > 0) {
            lapply(goRes, function(listOfHyperGs) {
                # list might contain results for MF, BP, CC. concatenate all
                p.value <-
                    do.call(c, lapply(listOfHyperGs, GOstats::pvalues))
                adj.p.value <- p.adjust(p.value, method = "BH")
                odds.ratio <-
                    do.call(c,
                            lapply(listOfHyperGs, GOstats::oddsRatios))
                expected.genes <-
                    do.call(c,
                            lapply(listOfHyperGs,
                                   GOstats::expectedCounts))
                matched.bicluster.ids <-
                    do.call(c,
                            lapply(listOfHyperGs,
                                   Category::geneIdsByCategory))
                matched.dataset.ids <-
                    do.call(c,
                            lapply(listOfHyperGs,
                                   Category::geneIdUniverse))
                matched.bicluster.genes <-
                    unlist(lapply(matched.bicluster.ids,
                                  length))
                matched.dataset.genes <-
                    unlist(lapply(matched.dataset.ids, length))
                names <-
                    suppressMessages(select(
                        GO.db,
                        keys = names(p.value),
                        columns = c("TERM"),
                        keytype = "GOID"
                    ))

                df <-
                    data.frame(
                        Name = names[, 2],
                        matched.bicluster.genes,
                        matched.dataset.genes,
                        expected.genes,
                        odds.ratio,
                        p.value,
                        adj.p.value
                    )
                df$matched.bicluster.ids <-
                    matched.bicluster.ids
                df$matched.dataset.ids <- matched.dataset.ids
                df
            })
        } else
            list())
    })

    goSummary <- reactive({
        goRes <- values$goRes
        if (length(goRes) > 1) {
            # summary data is same regardless of bicluster tested
            listOfHyperGs <- goRes[[1]]
            # list of matched IDs for each GO term
            matchedDatasetIds <- do.call(c,
                                         lapply(listOfHyperGs,
                                                Category::geneIdUniverse))
            # how many in the dataset were actually GO-annotated
            matchedDatasetSize <-
                length(unique(unlist(matchedDatasetIds)))
            # get the lowest adj. p-value for each bicluster
            significance <-
                vapply(
                    goDF(),
                    FUN.VALUE = numeric(1),
                    FUN = function(df)
                        min(df$adj.p.value)
                )
            if (any(significance < .Machine$double.eps)) {
                significance <- significance + .Machine$double.eps
            }
            significance <-
                sort(-log10(significance), decreasing = TRUE)

            return(list(
                universeSize = matchedDatasetSize,
                significance = significance
            ))
        }
    })

    # Plot of bicluster significance values
    output$goSigPlot <- renderPlot({
        goRes <- values$goRes
        validate(need(
            length(goRes) > 0 && !is.null(names(goRes)),
            "Please test for functional enrichment"
        ))

        bp <-
            barplot(
                height = goSummary()$significance,
                xaxs = "i",
                yaxs = "i",
                xaxt = "n",
                ylab = "-log10[Adj. p-value]"
            )
        las <-
            if (any(nchar(names(goRes)) > 5))
                2
        else
            1 # labels rotated?
        axis(
            side = 1,
            at = bp,
            pos = 0,
            labels = names(goSummary()$significance),
            las = las
        )
        abline(h = -log10(c(0.05)),
               col = "#2ca25f",
               lty = "dashed")
    })

    # Get input bicluster to display results for
    output$goBicluster <- renderUI({
        choices <- if (length(goDF()) > 0)
            names(goDF())
        else
            list()
        val <- if (length(choices) > 0)
            choices[[1]]
        else
            NULL
        selectInput(
            "goBicluster",
            label = "Select bicluster:",
            choices = choices,
            selected = val,
            multiple = FALSE
        )
    })

    # Display detailed GO results for one bicluster
    output$goTermTable <- DT::renderDT({
        validate(need(
            length(values$goRes) > 0 && !is.null(names(values$goRes)) &&
                !is.null(input$goBicluster),
            "Please test GO enrichment"
        ))
        df <-
            goDF()[[input$goBicluster]][, 1:7] # don't include gene lists
        return(
            DT::datatable(
                df,
                options = list(paging = FALSE, info = FALSE),
                colnames = c(
                    "GO ID",
                    "Name",
                    "Matched bicluster genes",
                    "Matched dataset genes",
                    "Expected genes",
                    "Odds ratio",
                    "p-value",
                    "Adj. p-value"
                ),
                fillContainer = TRUE,
                autoHideNavigation = TRUE,
                selection = 'single'
            )
        )
    })

    observeEvent(
        input$goTabGenes,
        {selection <- input$goTermTable_rows_selected
        if (length(selection) > 0) {
            # assumes goDF() has not changed since goTermTable was rendered
            selection <-
                row.names(goDF()[[input$goBicluster]])[selection]

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
        if (!is.null(bicluster) && bicluster %in% names(goDF())) {
            choices <- row.names(goDF()[[bicluster]])
            updateSelectInput(
                session,
                inputId = "goTerm",
                choices = choices,
                selected = choices[1]
            )
        }
    })

    # FIXME these get too large; can we put them in a scrolling div?
    output$goBiclusterGenes <- renderText({
        validate(need(length(goDF()) > 0, "Please test GO enrichment"))
        validate(need(!is.null(input$goTerm), "Please select a goTerm"))
        validate(need(
            input$goBicluster %in% names(goDF()),
            "Please select a valid bicluster"
        ))

        genes <-
            goDF()[[input$goBicluster]][input$goTerm, "matched.bicluster.ids"]
        paste(unlist(genes), collapse = "\n")
    })

    output$goUniverseGenes <- renderText({
        validate(need(length(goDF()) > 0, "Please test GO enrichment"))
        validate(need(!is.null(input$goTerm), "Please select a goTerm"))
        validate(need(
            input$goBicluster %in% names(goDF()),
            "Please select a valid bicluster"
        ))

        genes <-
            goDF()[[input$goBicluster]][input$goTerm, "matched.dataset.ids"]
        paste(unlist(genes), collapse = "\n")
    })

    #### SESSION ####
    # stop App on closing the browser
    session$onSessionEnded(function() {
        stopApp()
    })

}
