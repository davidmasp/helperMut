
<!-- README.md is generated from README.Rmd. Please edit that file -->

# helperMut <a href=''><img src='helperMut.png' align="right" height="165" /></a>

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/davidmasp/helperMut.svg?branch=develop)](https://travis-ci.com/davidmasp/helperMut)
[![codecov](https://codecov.io/gh/davidmasp/helperMut/branch/develop/graph/badge.svg?token=9jkMksb2mk)](https://codecov.io/gh/davidmasp/helperMut)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

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
#> helperMut Coverage: 34.92%
#> R/plots.R: 0.00%
#> R/utils.R: 2.86%
#> R/profiles.R: 12.82%
#> R/genome.R: 16.47%
#> R/indels.R: 22.22%
#> R/regions.R: 41.51%
#> R/muts.R: 61.27%
```

## Usage

### Simplify mutation codes

``` r
mutations = c("TCA>T","TGA>A")
simplify_muts(mutations)
#> [1] "TCA>T" "TCA>T"
```

### Cosine similarities from Cosmic signatures

``` r
# generate 3 random sigs
set.seed(42)
de_novo_sig = runif(n = 96*3,min = 0,max = 1)
de_novo_sig = matrix(de_novo_sig, ncol = 3)
rownames(de_novo_sig) = pos_ms96
colnames(de_novo_sig) = LETTERS[1:3]

cosmic_sigs = download_signature_set(type = "cosmic")

res = compare_signature_sets(x = de_novo_sig,y = cosmic_sigs)
heatmap(res)
```

![](README-unnamed-chunk-3-1.png)<!-- -->

### Colors and Plots

The standard colors used for mutational profiles in papers are available
as a variable.

``` r
tr_colors
#>       C>G       C>A       C>T       A>T       A>G       A>C 
#> "#000000" "#00ceff" "#ff2926" "#c9c9c9" "#95d230" "#ffbebe"

tr_colors_ct = tr_colors
names(tr_colors_ct) = simplify_muts(names(tr_colors_ct),
                                    simplify_set = c("C","T"))

tr_colors_ct
#>       C>G       C>A       C>T       T>A       T>C       T>G 
#> "#000000" "#00ceff" "#ff2926" "#c9c9c9" "#95d230" "#ffbebe"
```

## Other relevant packages

Some functionalities implemented in this packages are also available in
other packages from bioconductor. Hopefully, the current implementation
in this package will have an added value to the user.

  - [SomaticSignatures](http://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html)
    It extracts Mutation Subtype information from VCFs (restricted to a
    defined ctx). It also performs NMF/PCA signature extraction which is
    not covered in helperMut.
  - [MutationalPatterns](http://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html)
    Same as somaticSignatures.

## Contribute

If you want to contribute to the package, please fork the repo and
submit a PR. Currently the package is under development so no features
are explicitily requested.
