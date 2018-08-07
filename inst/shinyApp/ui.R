tagList(
  shinyjs::useShinyjs(),
  tags$head(tags$style(type="text/css", ".tab-content { height: 82vh;
                       !important; }" ),
            tags$script('
                              var dimension = [0, 0];
                              $(document).on("shiny:connected", function(e) {
                              dimension[0] = window.innerWidth;
                              dimension[1] = window.innerHeight;
                              Shiny.onInputChange("dimension", dimension);
                              });
                              $(window).resize(function(e) {
                              dimension[0] = window.innerWidth;
                              dimension[1] = window.innerHeight;
                              Shiny.onInputChange("dimension", dimension);
                              });
                              ')),
  #### UI ##################################################################
  navbarPage(
    theme = shinythemes::shinytheme("yeti"), inverse = TRUE, "mfBiclust UI",
    #### Data tabpanel ####
    source(system.file('shinyApp', "fileUI.R", package='mfBiclust'), local = TRUE)$value,
    #### Bicluster tabpanel ####
    source(system.file('shinyApp', "biclusterUI.R", package='mfBiclust'), local = TRUE)$value,
    source(system.file('shinyApp', "bcvUI.R", package='mfBiclust'), local = TRUE)$value, 
    source(system.file('shinyApp', "goUI.R", package='mfBiclust'), local = TRUE)$value,
    id = "navbarpage", fluid = TRUE,
    footer = tags$footer(verbatimTextOutput("status"))
  )
)