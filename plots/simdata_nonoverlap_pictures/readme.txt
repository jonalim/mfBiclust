BCV[ESA] 3 iterations ESA. Per authors instructions, k* = argmin[k] summation( PE, 100 randomized holdouts ) / 100
BCV[PCA] k* = argmin[k] summation( ||residuals||, 10x10 fold holdouts )

Both algorithms run 40 replications, both median and mean k* were examined. Median was slightly more accurate for BCV[PCA].