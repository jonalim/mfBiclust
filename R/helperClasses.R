# inner object
setClass("genericFactorization", slots = list(W = "matrix", H = "matrix"))

# outer object
setClass("genericFit", slots = list(fit = "genericFactorization", 
                                    method = "character"))