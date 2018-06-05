setClass("genericFactorization", slots = list(W = "matrix", H = "matrix"))

setClass("genericFit", slots = list(fit = "genericFactorization", 
                                    method = "character"))