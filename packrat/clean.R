packrat::set_opts("auto.snapshot" = FALSE)
up <- packrat::unused_packages()
up <- up[sapply(up, function(x) ! x$name %in% c("BH", "plogr", "Rcpp", "BiocGenerics"))]
packrat::clean(up)

if(!requireNamespace("devtools")) {
  install.packages("devtools")
}
if(packageVersion("devtools") < package_version("1.13.5.9000")) {
  devtools::install_github("r-lib/devtools", ref = "33a9404583833f162f84b6a931cca2ccca2dd8d7")
}

devtools::install_version("roxygen2", version = "6.0.1")
