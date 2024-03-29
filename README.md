MR.CUE
=======
  
  **MR.CUE** is a package for Mendelian randomization informs shared genetic
etiology underlying exposure and outcome by
interrogating correlated horizontal pleiotropy.

Installation
============
  Install the development version of **MR.CUE** by use of the 'devtools' package. Note that **MR.CUE** depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```
library(devtools)
install_github("QingCheng0218/MR.CUE@main")
```

If you have errors installing this package on Linux, try the following commands in R:
  ```
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252") 
Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

library(devtools)
install_github("QingCheng0218/MR.CUE@main")
```

Usage
=========
  The ['MR.CUE' vignette](https://github.com/QingCheng0218/MR.CUE/blob/main/vignettes/MR-CUE.pdf) will provide a good start point for Mendelian randomization analysis using **MR.CUE** package. 

References
==========

Qing Cheng, Xiao Zhang, Lin S. Chen<sup>+</sup>, Jin Liu<sup>+</sup>. Mendelian randomization accounting for complex correlated horizontal pleiotropy while elucidating shared genetic etiology. Nat Commun 13, 6490 (2022). https://doi.org/10.1038/s41467-022-34164-1

Development
===========
  
  This package is developed and maintained by Qing Cheng (qingcheng0218@gmail.com). 
  Cite the code:[![DOI](https://zenodo.org/badge/509940910.svg)](https://zenodo.org/badge/latestdoi/509940910)
