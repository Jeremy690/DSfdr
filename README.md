# DSfdr
## Overview
This R package selects variables in linear regression and Gaussian graphical model with FDR under control.

## Running examples
- You can find package **vignette** in the vignettes fold


## Installation

You can install the package using 

```R
install.packates('devtools')
devtools::install_github("Jeremy690/DSfdr/DSfdr")
```

To install with vignettes, using 

```R
devtools::install_github("Jeremy690/DSfdr/DSfdr", build_vignettes = TRUE)
```

The vignettes takes approximately 5 mins to knit.


## Feedback

If you encounter error or would like to provide feedback, please use [Github -> Issues](https://github.com/LinBuyu/DSfdr/issues) to reach us. Thank you! 


## Reproducible File
We have uploaded an Rmd file to reproduce all the results in the paper [False Discovery Rate Control via Data Splitting](https://arxiv.org/pdf/2002.08542.pdf). The knockoff. R and MBHq.R are used to compare our methods with other two methods for linear regression. The data.txt are used in one of our simulation results. 
