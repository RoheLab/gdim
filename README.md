


<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `gdim`: Graph dimension

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`gdim` is an R package for estimating graph dimension `k` (e.g., the
number of communities in a network). The key function `eigcv` uses the
cross-validated eigenvalues to test the statistical significance of
sample eigenvectors. The test statistic enjoys a simple central limit
theorem. As such, `eigcv` provides p-values for individual hypothesized
`k`.

## Installation

You can install `gdim` from [GitHub](https://github.com/) with:

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
A <- sample_sparse(sbm, poisson_edges = F, allow_self_loops = F)
```

Here, `A` is the adjacency matrix.

Now, we call the `eigcv` function with `k_max=10` to estimate graph
dimension.

``` r
eigcv(A, k_max = 10)
#> Estimated graph dimension:    5
#> 
#> Number of bootstraps:         10
#> Edge splitting probabaility:  0.1
#> Significance level:       0.05
#> 
#>  ------------ Summary of Tests ------------
#>   k       zbar      zmin         pval         padj
#>   1 62.2747913 61.589009 0.000000e+00 0.000000e+00
#>   2 10.9163348  9.451899 4.815205e-28 4.815205e-28
#>   3 10.5640433  9.215151 2.187029e-26 2.187029e-26
#>   4  9.3309331  7.086002 5.247143e-21 5.247143e-21
#>   5  6.4823670  4.474127 4.514734e-11 4.514734e-11
#>   6 -1.2957003 -2.066368 9.024606e-01 9.024606e-01
#>   7 -0.7079412 -2.278473 7.605091e-01 7.605091e-01
#>   8 -1.1075418 -2.626536 8.659701e-01 8.659701e-01
#>   9 -1.3166394 -3.038380 9.060202e-01 9.060202e-01
#>  10  0.1547320 -2.244862 4.385163e-01 4.385163e-01
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
