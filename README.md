# PopGenBounds 

<!-- badges: start -->

[![R-CMD-check](https://github.com/nalcala/PopGenBounds/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nalcala/PopGenBounds/actions/workflows/R-CMD-check.yaml) [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) [![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/license/gpl-3-0-only)

<!-- badges: end -->

The *PopGenBounds* R package provides functions and visualization to explore diversity and differentiation statistics (FST, Hedrick's G'ST, and Jost's D) with respects to their mathematical bounds as functions of the frequency of the most frequent allele.

## Installation

You can install this package from GitHub with the following R code:

``` r
# install.packages("devtools")
devtools::install_github("nalcala/PopGenBounds")
```

## Usage

We load the package

``` r
library(PopGenBounds)
```

Let's first create an allele frequency matrix with subpopulations as rows and alleles as columns:

``` r
popa = matrix(c(0.4,0.4,0.2,rep(0,3*5),
                rep(0,3),0.4,0.4,0.2,rep(0,3*4),
                rep(0,3*2),0.4,0.4,0.2,rep(0,3*3),
                rep(0,3*3),0.4,0.4,0.2,rep(0,3*2),
                rep(0,3*4),0.4,0.4,0.2,rep(0,3),
                rep(0,3*5),0.4,0.4,0.2),ncol=18,byrow = T)
```

The Diff function computes the differentiation statistics along with their bounds given the most frequent allele *M*

```         
Diff(popa)
```

We can plot the allele frequency table with function ggfreqtable:

``` r
ggfreqtable(popa)
```

And plot the value within the bounds of FST, G'ST, and D with function ggbounds:

``` r
ggbounds(M=Diff(popa)$value[1],FST=Diff(popa)$value[2],GpST = Diff(popa)$value[3],D=Diff(popa)$value[4],K=nrow(popa))
```

## Citations

Alcala, N., Lim, Y., & Rosenberg, N. A. (In prep)

Alcala, N., & Rosenberg, N. A. (2022). Mathematical constraints on F ST: multiallelic markers in arbitrarily many populations. Philosophical Transactions of the Royal Society B, 377(1852), 20200414.

Alcala, N., & Rosenberg, N. A. (2019). G'ST, Jost's D, and FST are similarly constrained by allele frequencies: A mathematical, simulation, and empirical study. Molecular ecology, 28(7), 1624-1636.

Alcala, N., & Rosenberg, N. A. (2017). Mathematical constraints on F ST: biallelic markers in arbitrarily many populations. Genetics, 206(3), 1581-1600.
