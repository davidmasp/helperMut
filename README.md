
<!-- README.md is generated from README.Rmd. Please edit that file -->

# helperMut <a href=''><img src='helperMut.png' align="right" height="165" /></a>

[![codecov](https://codecov.io/gh/davidmasp/helperMut/branch/develop/graph/badge.svg?token=9jkMksb2mk)](https://codecov.io/gh/davidmasp/helperMut)

## Overview

HyperMut provides basic functionalities and helper functions to work
with mutational and genomic data in R. It builds up on pre-existing
tools from bioc to set up a framework for analysis. I mostly use it for
my work in cancer genomics but in principle any other field working in
mutations or population diversity could be used.

## Installation

``` r
## If still private, make sure that there is a enviroment variable (.Renviron)
## with the token (PAT). The name of the variable should be GITHUB_PAT
remotes::install_github("davidmasp/helpermut@develop")
```

``` r

covr::package_coverage()
#> helperMut Coverage: 48.42%
#> R/plots.R: 0.00%
#> R/utils.R: 2.86%
#> R/genome.R: 16.47%
#> R/regions.R: 41.51%
#> R/indels.R: 76.47%
#> R/muts.R: 78.32%
```

## Usage

### Simplify mutation codes

``` r
mutations = c("TCA>T","TGA>A")
simplify_muts(mutations)
#> [1] "TCA>T" "TCA>T"
```

## Other relevant packages

## Contribute
