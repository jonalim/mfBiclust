bce <- BiclusterExperiment(yeast_benchmark[[1]])
addStrat(bce, k = 3)
\donttest{biclusterGUI(bce)}
