dir = "../biclustlib-master/biclustlib/datasets/data/yeast_benchmark/"
files <- list.files(path = dir)
yeast_benchmark <- lapply(files, function(x) {
  tab <- read.table(file = paste0(dir, x), sep = "",
                    header = FALSE, row.names = 2, stringsAsFactors = FALSE)
  tab[, 2:ncol(tab)]
})
usethis::use_data(yeast_benchmark)

dir = "../biclustlib-master/biclustlib/datasets/data/cancer_benchmark/"
files <- c("chen-2002_database.txt",
           "chowdary-2006_database.txt",
           "nutt-2003-v3_database.txt",
           "pomeroy-2002-v1_database.txt",
           "nutt-2003-v2_database.txt",
           "singh-2002_database.txt",
           "alizadeh-2000-v1_database.txt",
           "dyrskjot-2003_database.txt",
           "liang-2005_database.txt",
           "pomeroy-2002-v2_database.txt",
           "west-2001_database.txt",
           "shipp-2002-v1_database.txt",
           "nutt-2003-v1_database.txt")
cancer_benchmark <- lapply(
  files, 
  function(x) {
    # Read without header because duplicate column names are not allowed
    tab <- read.table(file = paste0(dir, x), sep = "",
                      header = FALSE, stringsAsFactors = FALSE)
    # first row contains labels; convert to classification matrix
    labels <- factor(as.character(tab[1, 2:ncol(tab)]))
    classes <- levels(labels)
    labels <- do.call(rbind, lapply(labels, function(lab) {
      vec <- rep(0, length(classes))
      names(vec) <- classes
      vec[lab] <- 1
      vec
    }))
    
    # drop first row and column and convert all dataframe columns to numeric
    # do not include transcript identifiers (some transcript identifiers are
    # duplicates; skip dealing with those)
    tab <- sapply(tab[2:nrow(tab), 2:ncol(tab)], as.numeric)

    # for each dataset, include expression dataframe and classification matrix
    list(data = tab, labels = labels)
  }
)
usethis::use_data(cancer_benchmark, overwrite = TRUE, compress = "bzip2")
tools::checkRdaFiles("data/cancer_benchmark.rda")
