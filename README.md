


<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `gdim`: Graph dimensionality

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`gdim` is an R package for estimating graph dimensionality `k` (e.g.,
the number of communities in a network). The key function `eigcv` uses
the cross-validated eigenvalues to test the statistical significance of
sample eigenvectors. The test statistics enjoy a simple central limit
theorem. As such, `eigcv` provides p-values for individual hypothesized
`k`.

## Installation

You can install from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("RoheLab/gdim")
```

<!-- You can install the released version of epca from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("gdim") -->
<!-- ``` -->
<!-- or the development version from [GitHub](https://github.com/) with: -->
<!-- ``` r -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("RoheLab/gdim") -->
<!-- ``` -->

## Example

`eigcv()` is the main function in `gdim`. The function is simple to use
and easy to configurable. The single required parameter for the function
is the maximum possible dimension, `k_max`.

In the following example, we generate a random graph from the stochastic
block model (SBM) with 1000 nodes and 5 blocks.

``` r
library(fastRG)
B = matrix(0.1, 5, 5)
diag(B) = 0.3
sbm <- sbm(n = 1000, k = 5, B = B, expected_degree = 40)
A <- sample_sparse(sbm, poisson_edges = T, allow_self_loops = F)
```

Here, `A` is the adjacency matrix.

Now, we call the `eigcv` function with `k_max=10` to estimate graph
dimensionality.

``` r
eigcv(A, k_max = 10)
#> Estimated graph dimensionality:   5
#> 
#> Number of bootstraps:         10
#> Edge splitting probabaility:  0.1
#> Significance level:       0.05
#> 
#>  ------------ Summary of Tests ------------
#>   k       zbar      zmin         pval         padj
#>   1 62.0801262 61.515047 0.000000e+00 0.000000e+00
#>   2 10.0745208  8.086710 3.580449e-24 3.580449e-24
#>   3  9.9983025  8.410203 7.751587e-24 7.751587e-24
#>   4  8.5536034  6.424442 5.965218e-18 5.965218e-18
#>   5  9.4837958  8.118961 1.225988e-21 1.225988e-21
#>   6 -0.5617609 -1.763343 7.128605e-01 7.128605e-01
#>   7 -1.5824039 -3.940681 9.432213e-01 9.432213e-01
#>   8 -0.6574413 -4.045873 7.445514e-01 7.445514e-01
#>   9 -1.3668588 -2.862552 9.141652e-01 9.141652e-01
#>  10 -1.3150209 -2.685515 9.057486e-01 9.057486e-01
```

In this example, `eigcv()` suggests to choose `k=5`.

<!-- For more examples, please see the vignette:  -->
<!-- ```{r, eval=FALSE} -->
<!-- vignette("epca") -->
<!-- ``` -->

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/RoheLab/gdim/issues).

## Reference

Chen F, Roch S, Rohe K, & Yu S. “Estimating graph dimension with
cross-validated eigenvalues.” *In preparation*.
