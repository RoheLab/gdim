


<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `gdim`: Graph dimension

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`gdim` is an R package for estimating graph dimension `k` (e.g., the
number of communities in a network). The key function `eigcv` uses the
cross-validated eigenvalues to test the statistical significance of
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
#>   1 62.2108768 61.671218 0.000000e+00 0.000000e+00
#>   2 10.2355302  8.716879 6.872015e-25 6.872015e-25
#>   3 10.9104976  9.085051 5.134636e-28 5.134636e-28
#>   4  9.6183780  6.874632 3.343868e-22 3.343868e-22
#>   5  8.9396277  7.405616 1.952418e-19 1.952418e-19
#>   6 -0.2306637 -2.036427 5.912120e-01 5.912120e-01
#>   7 -1.6699556 -3.609351 9.525359e-01 9.525359e-01
#>   8 -1.4463405 -4.013657 9.259591e-01 9.259591e-01
#>   9 -0.6664051 -2.576775 7.474239e-01 7.474239e-01
#>  10 -0.7495700 -2.479426 7.732431e-01 7.732431e-01
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
