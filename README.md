


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
block model (SBM) with 1000 nodes and 5 blocks (as such, we would expect
the estimated graph dimension to be 5).

``` r
library(fastRG)
B = matrix(0.1, 5, 5)
diag(B) = 0.3
sbm <- sbm(n = 1000, k = 5, B = B, expected_degree = 40)
A <- sample_sparse(sbm, poisson_edges = F, allow_self_loops = F)
```

Here, `A` is the adjacency matrix.

Now, we call the `eigcv()` function with `k_max=10` to estimate graph
dimension.

``` r
eigcv_result <- eigcv(A, k_max = 10)
eigcv_result
#> Estimated graph dimension:    5
#> 
#> Number of bootstraps:         10
#> Edge splitting probabaility:  0.1
#> Significance level:       0.05
#> 
#>  ------------ Summary of Tests ------------
#>   k          z        pvals         padj
#>   1 41.7764210 1.000000e-32 1.000000e-32
#>   2  7.3449938 1.028843e-13 1.028843e-13
#>   3  6.3361165 1.178143e-10 1.178143e-10
#>   4  5.8773752 2.084113e-09 2.084113e-09
#>   5  6.1895117 3.017544e-10 3.017544e-10
#>   6 -1.0575585 8.548716e-01 8.548716e-01
#>   7 -1.0084801 8.433880e-01 8.433880e-01
#>   8 -0.9298564 8.237773e-01 8.237773e-01
#>   9 -0.8521407 8.029320e-01 8.029320e-01
#>  10 -0.9469157 8.281591e-01 8.281591e-01
```

In this example, `eigcv()` suggests to choose `k=5`.

<!-- For more examples, please see the vignette:  -->
<!-- ```{r, eval=FALSE} -->
<!-- vignette("gdim") -->
<!-- ``` -->

To visualize the result, use `plot()` which returns a `ggplot` object.
The function displays the test statistic (z score) for each hypothesized
graph dimension. For more options, check out the manual.

<!-- ```{r plot_eigcv} -->
<!-- plot(eigcv_result) -->
<!-- ``` -->

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on our GitHub [Issuas
page](https://github.com/RoheLab/gdim/issues).

## Reference

Chen F, Roch S, Rohe K, & Yu S. “Estimating graph dimension with
cross-validated eigenvalues.” *In preparation*.
