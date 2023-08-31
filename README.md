
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gdim

<!-- badges: start -->
<!-- badges: end -->

`gdim` estimates graph dimension using cross-validated eigenvalues, via
the graph-splitting technique developed in
<http://arxiv.org/abs/2108.03336>. Theoretically, the method works by
computing a special type of cross-validated eigenvalue which follows a
simple central limit theorem. This allows users to perform hypothesis
tests on the rank of the graph.

## Installation

You can install `gdim` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("RoheLab/gdim")
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
eigcv_result
#> Estimated graph dimension:    5
#> 
#> Number of bootstraps:         10
#> Edge splitting probabaility:  0.1
#> Significance level:       0.05
#> 
#>  ------------ Summary of Tests ------------
#>   k          z        pvals         padj
#>   1 60.1876947 1.000000e-32 1.000000e-32
#>   2 12.3026816 1.000000e-32 1.000000e-32
#>   3 11.7647552 2.964906e-32 2.964906e-32
#>   4 11.6845190 7.647145e-32 7.647145e-32
#>   5  9.3815651 3.250029e-21 3.250029e-21
#>   6 -1.1955845 8.840706e-01 8.840706e-01
#>   7 -1.9103184 9.719539e-01 9.719539e-01
#>   8 -1.8069433 9.646144e-01 9.646144e-01
#>   9 -1.2305367 8.907519e-01 8.907519e-01
#>  10 -0.3782357 6.473723e-01 6.473723e-01
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
