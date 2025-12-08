# Print cross-validated eigenvalues

Print cross-validated eigenvalues

## Usage

``` r
# S3 method for class 'eigcv'
print(x, ...)
```

## Arguments

- x:

  An `eigcv` object created by a call to
  [`eigcv()`](https://rohelab.github.io/gdim/dev/reference/eigcv.md).

- ...:

  Ignored.

## Value

`x`, but invisibly.

## Examples

``` r
library(fastRG)

set.seed(27)

B <- matrix(0.1, 5, 5)
diag(B) <- 0.3

model <- sbm(
  n = 200,
  B = B,
  expected_degree = 20,
  poisson_edges = FALSE,
  allow_self_loops = FALSE
)

A <- sample_sparse(model)

eigs<- eigcv(A, k_max = 7)
eigs
#> Estimated graph dimension:    1
#> 
#> Number of bootstraps:         10
#> Edge splitting probabaility:  0.1
#> Significance level:       0.05
#> 
#>  ------------ Summary of Tests ------------
#>  k           z        pvals         padj
#>  1 11.85408886 1.024706e-32 1.024706e-32
#>  2  0.17249449 4.315244e-01 4.315244e-01
#>  3  0.54647185 2.923708e-01 2.923708e-01
#>  4 -0.85440917 8.035608e-01 8.035608e-01
#>  5 -0.09633109 5.383712e-01 5.383712e-01
#>  6  0.09554781 4.619399e-01 4.619399e-01
#>  7 -0.29814636 6.172043e-01 6.172043e-01
#> 

plot(eigs, type = "z-score")    # default

plot(eigs, type = "adjacency")

plot(eigs, type = "laplacian")


```
