# gdim: Estimate Graph Dimension using Cross-Validated Eigenvalues

Cross-validated eigenvalues are estimated by splitting a graph into two
parts, the training and the test graph. The training graph is used to
estimate eigenvectors, and the test graph is used to evaluate the
correlation between the training eigenvectors and the eigenvectors of
the test graph. The correlations follow a simple central limit theorem
that can be used to estimate graph dimension via hypothesis testing, see
Chen et al. (2021) [arXiv:2108.03336](https://arxiv.org/abs/2108.03336)
for details.

## See also

Useful links:

- <https://github.com/RoheLab/gdim>

- <https://rohelab.github.io/gdim/>

- Report bugs at <https://github.com/RoheLab/gdim/issues>

## Author

**Maintainer**: Alex Hayes <alexpghayes@gmail.com>
([ORCID](https://orcid.org/0000-0002-4985-5160)) \[copyright holder\]

Authors:

- Fan Chen <fchen365@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-4508-6023))

- Karl Rohe <karlrohe@wisc.edu>
