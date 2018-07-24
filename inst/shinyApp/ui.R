tagList(
  shinyjs::useShinyjs(),
  #### UI ##################################################################
  navbarPage(
    theme = shinythemes::shinytheme("yeti"), inverse = TRUE, "mfBiclust UI",
    #### Data tabpanel ####
    source(system.file('shinyApp', "fileUI.R", package='mfBiclust'), local = TRUE)$value,
    #### Bicluster tabpanel ####
    source(system.file('shinyApp', "biclusterUI.R", package='mfBiclust'), local = TRUE)$value,
    source(system.file('shinyApp', "bcvUI.R", package='mfBiclust'), local = TRUE)$value, 
    source(system.file('shinyApp', "goUI.R", package='mfBiclust'), local = TRUE)$value,
    id = "navbarpage", fluid = TRUE)
)