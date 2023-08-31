
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gdim

<!-- badges: start -->

[![R-CMD-check](https://github.com/RoheLab/gdim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RoheLab/gdim/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/RoheLab/gdim/branch/main/graph/badge.svg)](https://app.codecov.io/gh/RoheLab/gdim?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/gdim)](https://CRAN.R-project.org/package=gdim)
<!-- badges: end -->

`gdim` estimates graph dimension using cross-validated eigenvalues, via
the graph-splitting technique developed in
<http://arxiv.org/abs/2108.03336>. Theoretically, the method works by
computing a special type of cross-validated eigenvalue which follows a
simple central limit theorem. This allows users to perform hypothesis
tests on the rank of the graph.

## Installation

You can install `gdim` from CRAN with:

``` r
install.packages("gdim")

# to get the development version from GitHub:
install.packages("pak")
pak::pak("RoheLab/gdim")
```

## Example

`eigcv()` is the main function in `gdim`. The single required parameter
for the function is the maximum possible dimension, `k_max`.

In the following example, we generate a random graph from the stochastic
block model (SBM) with 1000 nodes and 5 blocks (as such, we would expect
the estimated graph dimension to be 5).

``` r
library(fastRG)
#> Loading required package: Matrix

B <- matrix(0.1, 5, 5)
diag(B) <- 0.3

model <- sbm(
  n = 1000,
  k = 5,
  B = B,
  expected_degree = 40,
  poisson_edges = FALSE,
  allow_self_loops = FALSE
)

A <- sample_sparse(model)
```

Here, `A` is the adjacency matrix.

Now, we call the `eigcv()` function with `k_max=10` to estimate graph
dimension.

``` r
library(gdim)

eigcv_result <- eigcv(A, k_max = 10)
#> 'as(<dsCMatrix>, "dgCMatrix")' is deprecated.
#> Use 'as(., "generalMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
eigcv_result
#> Estimated graph dimension:    5
#> 
#> Number of bootstraps:         10
#> Edge splitting probabaility:  0.1
#> Significance level:       0.05
#> 
#>  ------------ Summary of Tests ------------
#>   k          z        pvals         padj
#>   1 60.1118360 2.220446e-16 2.220446e-15
#>   2 12.4766591 2.220446e-16 2.220446e-15
#>   3 12.0401480 2.220446e-16 2.220446e-15
#>   4 10.7393316 2.220446e-16 2.220446e-15
#>   5  8.9024334 2.220446e-16 2.220446e-15
#>   6 -1.1480263 8.745212e-01 1.000000e+00
#>   7 -1.6158410 9.469357e-01 1.000000e+00
#>   8 -0.9936334 8.397993e-01 1.000000e+00
#>   9 -0.8575451 8.044281e-01 1.000000e+00
#>  10 -1.1543480 8.758212e-01 1.000000e+00
```

In this example, `eigcv()` suggests `k=5`.

To visualize the result, use `plot()` which returns a `ggplot` object.
The function displays the test statistic (z score) for each hypothesized
graph dimension.

``` r
plot(eigcv_result)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## Reference

Chen, Fan, Sebastien Roch, Karl Rohe, and Shuqi Yu. “Estimating Graph
Dimension with Cross-Validated Eigenvalues.” ArXiv:2108.03336 \[Cs,
Math, Stat\], August 6, 2021. <http://arxiv.org/abs/2108.03336>.
