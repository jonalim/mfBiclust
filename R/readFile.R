#' @export
chooseFile <- function() {
  input <- NULL
  while(is.null(input)) {
    
    filepath <- file.choose()
    message("Parsing your file...")
    expanded <- strsplit(filepath, '\\.')[[1]]
    ext <- expanded[length(expanded)]
    if(ext == "xls" || ext == "xlsx") {
      tryCatch({
        input <- readxl::read_excel(filepath, 
                                    sheet = "Sheet1", col_names = FALSE)
      }, warning = function(w) {
        res <- tcltk::tkmessageBox(message = paste(w, "\nPlease click OK to proceed"), icon = "warning", 
                                   type = "okcancel", default = "ok")
        if(as.character(res) == "cancel") {input == NULL}
        else {input <- readxl::read_excel(filepath, 
                                          sheet = "Sheet1", col_names = FALSE)}
      }, error = function(e) {
        res <- tcltk::tkmessageBox(message = paste(e, "\nmfBiclust cannot continue"), icon = "error", type = "retrycancel", default = "retry")
        if(as.character(res) == "retry") {input == NULL} # All installed packages will go into this directory
        else {
          # q()
        }
      }
      )
    }
    else if(ext == "csv") {
      tryCatch({
        # CSV auto-header detection isn't working? Got 29 rows when 30 expected
        input <- read.table(filepath, sep = ",", quote = "\"",
                            fill = TRUE, comment.char = "", colClasses = "numeric")
      }, warning = function(w) {
        res <- tcltk::tkmessageBox(message = paste(w, "\nPlease click OK to proceed"), icon = "warning", 
                                   type = "okcancel", default = "ok")
        if(as.character(res) == "cancel") {input == NULL}
        else {input <- read.csv(filepath, colClasses = "numeric")}
      }, error = function(e) {
        res <- tcltk::tkmessageBox(message = paste(e, "\nmfBiclust cannot continue"), icon = "error", type = "retrycancel", default = "retry")
        if(as.character(res) == "retry") {input == NULL} # All installed packages will go into this directory
        else {
          # q()
        }
      }
      )
    } else {
      tryCatch({
        input <- read.table(filepath, colClasses = "numeric")
      }, warning = function(w) {
        res <- tcltk::tkmessageBox(message = paste(w, "\nPlease click OK to proceed"), icon = "warning", 
                                   type = "okcancel", default = "ok")
        if(as.character(res) == "cancel") {input == NULL}
        else {input <- read.table(filepath, colClasses = "numeric")}
      }, error = function(e) {
        res <- tcltk::tkmessageBox(message = paste(e, "\nmfBiclust cannot continue"), icon = "error", type = "retrycancel", default = "retry")
        if(as.character(res) == "retry") {input == NULL} # All installed packages will go into this directory
        else {
          # q()
        }
      }
      )
      
    }
    
  }
  
  as.matrix(input)
}
