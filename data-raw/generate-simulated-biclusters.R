invisible(duplicable(str = "genera"))
simdata_3_nonoverlap <- genSimData(n = 3, overlapped = FALSE, noise = 0.01)
simdata_3_overlap <- genSimData(n = 3, overlapped = TRUE, noise = 0.01)
simdata_10_nonoverlap <- genSimData(n = 10, overlapped = FALSE, noise = 0.01)
simdata_10 <- genSimData(n = 10, overlapped = TRUE, noise = 0.01)

simdata <- list(simdata_3_nonoverlap, simdata_3_overlap, simdata_10_nonoverlap,
                simdata_10)
<<<<<<< HEAD
names(simdata) <- c("simdata_3_nonoverlap", "simdata_3_overlap",
                    "simdata_10_nonoverlap", "simdata_10")
=======
names(simdata) <- c("nonoverlap3", "overlap3",
                    "nonoverlap10", "overlap10")
>>>>>>> Added simulation data and  to the correct folders
usethis::use_data(simdata, overwrite = TRUE)
