set.seed(12345)
nonoverlap_3biclusters <- genSimData(n = 3, overlapped = FALSE, noise = 0.01,
                                     dynamSize = FALSE, dimx = 40)
nonoverlap_15biclusters <- genSimData(n = 15, overlapped = FALSE, noise = 0.01,
                                      dynamSize = FALSE, dimx = 40)

nonoverlap_5height <- genSimData(n= 3, overlapped = FALSE, noise = 0.01,
                               dynamSize = FALSE, dimx = 200, dimy = 40)

nonoverlap_5width <- genSimData(n= 3, overlapped = FALSE, noise = 0.01,
                                 dynamSize = FALSE, dimx = 40, dimy = 200)

overlap_3biclusters <- genSimData(n = 3, overlapped = TRUE, noise = 0.01,
                                  dynamSize = FALSE, dimx = 40)
set.seed(123)
overlap_4biclusters <- genSimData(n = 4, overlapped = TRUE, noise = 0.01,
                                  dynamSize = FALSE, dimx = 40)

simdata <- list(nonoverlap_3biclusters, nonoverlap_15biclusters, nonoverlap_5height,
                nonoverlap_5width, overlap_3biclusters, overlap_4biclusters)
names(simdata) <- c("nonoverlap_3biclusters", "nonoverlap_15biclusters",
                    "nonoverlap_5height", "nonoverlap_5width",
                    "overlap_3biclusters", "overlap_4biclusters")
usethis::use_data(simdata, overwrite = TRUE)
