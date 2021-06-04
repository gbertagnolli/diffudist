
<!-- README.md is generated from README.Rmd. Please edit that file -->

# diffudist <img src="man/figures/diffudist.png" align="right" alt="" width="120"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/gbertagnolli/diffudist/workflows/R-CMD-check/badge.svg)](https://github.com/gbertagnolli/diffudist/actions)
<!-- badges: end -->

## Overview

The `diffudist` package provides several functions for evaluating the
diffusion distance between nodes of a complex network.

## Installation

``` r
# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("gbertagnolli/diffudist")
```

## Usage

Additionally to `diffudist` you will also need the `igraph` package,
because the main arguments of the functions in `diffudist` are networks
as `igraph` objects.

``` r
library(diffudist)
library(igraph)
```

### Examples

``` r
g <- sample_pa(n = 100, power = 1, m = 2)
D <- get_DDM(g, tau = 2, type = "Normalized Laplacian", verbose = FALSE)
#> Unweighted network.
#> Evaluating the Normalized Laplacian matrix
MERW_Pt <- get_diffu_Pt(g, tau = 2, type = )
#> Unweighted network.
#> Evaluating the Normalized Laplacian matrix
```

The probability transition matrix returned from
`get_diffusion_probability_matrix` (or its shortened version
`get_diffu_Pt`) is the matrix *e*<sup> − *τ**L*<sub>rw</sub></sup>. The
diffusion dynamics is controlled by the specific Laplacian matrix
*L*<sub>rw</sub> = *I* − *T*<sub>rw</sub>, where *T*<sub>rw</sub> is the
jump matrix of the discrete-time random walk corresponding to our
continuous-time dynamics.

## References

Bertagnolli, G., & De Domenico, M. (2021). *Diffusion geometry of
multiplex and interdependent systems*. Physical Review E, 103(4),
042301. [DOI:
10.1103/PhysRevE.103.042301](https://doi.org/10.1103/PhysRevE.103.042301),
[arXiv: 2006.13032](https://arxiv.org/abs/2006.13032),
[my-website](https://gbertagnolli.github.io/publication/ml-diffusion/).
