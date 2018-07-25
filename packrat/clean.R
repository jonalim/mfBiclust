packrat::set_opts("auto.snapshot" = FALSE)
up <- packrat::unused_packages()
up <- up[unlist(sapply(up, function(x) ! x$name %in% c("BH", "plogr", "BiocGenerics")))]
packrat::clean(up, force = TRUE)

try(unloadNamespace("devtools"))
if(!requireNamespace("devtools")) {
  install.packages("devtools")
}
if(packageVersion("devtools") < package_version("1.13.5.9000")) {
  devtools::install_github("r-lib/devtools", ref = "33a9404583833f162f84b6a931cca2ccca2dd8d7")
}

devtools::install_github("klutometis/roxygen#760")
devtools::document(roclets=c('rd', 'collate', 'namespace', 'vignette'))
devtools::install()
packrat::snapshot(ignore.stale = TRUE, prompt = TRUE)
