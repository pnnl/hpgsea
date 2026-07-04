
- [fast.ssgsea](#fastssgsea)
  - [Overview](#overview)
  - [Installation](#installation)
    - [macOS](#macos)
    - [Windows](#windows)
    - [Linux](#linux)
    - [Install](#install)
  - [Usage](#usage)
    - [Simulate Data](#simulate-data)
    - [Runtime and Results](#runtime-and-results)
    - [Session Information](#session-information)
  - [Benchmarking](#benchmarking)
    - [fast-ssGSEA](#fast-ssgsea)
    - [FGSEA-simple](#fgsea-simple)
  - [References](#references)

# fast.ssgsea

<!-- badges: start -->

[![R-CMD-check](https://github.com/pnnl/fast.ssgsea/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pnnl/fast.ssgsea/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**NOTICE:** While this R package was based on the
[ssGSEA2.0](https://github.com/broadinstitute/ssGSEA2.0) repository,
neither perform single-sample Gene Set Enrichment Analysis (ssGSEA) as
originally described by Barbie, *et al.* ([Barbie et al.
2009](#ref-barbie-systematic-2009)). They are instead modifications of
pre-ranked GSEA that calculate the enrichment score (ES) differently and
support testing directional gene sets (details below). The package and
fast-ssGSEA name will be changed in the future.

## Overview

`fast.ssgsea` is an R package ([R Core Team 2026](#ref-R-core-team)) for
a highly optimized variant of pre-ranked Gene Set Enrichment Analysis
(GSEA) ([Subramanian et al. 2005](#ref-subramanian-gene-2005)). Unlike
standard GSEA, fast-ssGSEA is capable of testing gene sets where each
gene has an expected direction of change (up- or down-regulation;
indicated by appending a “;u” or “;d” to the end of every gene in a set)
from a prior experiment.

fast-ssGSEA is based on Post-Translational Modification Signature
Enrichment Analysis (PTM-SEA) ([Krug et al.
2019](#ref-krug-curated-2019)), and it borrows optimization techniques
from the simple implementation of Fast Gene Set Enrichment Analysis
(FGSEA-simple) ([Korotkevich et al. 2021](#ref-korotkevich-fast-2021)).

The primary function, `fast_ssgsea`, accepts a vector of signed
statistics with genes or other molecules as names. The values must be
approximately symmetric around zero, with more extreme values indicating
greater importance. A named list of gene sets (more generally, molecular
signatures) is also required. Other arguments control the behavior of
fast-ssGSEA, and they are described in the function documentation.

The package also contains a `read_gmt` function, which reads a Gene
Matrix Transposed (GMT) file to construct a named list of gene sets for
use with `fast_ssgsea`.

## Installation

R version 4.0.0 or greater is required to install `fast.ssgsea`.

### macOS

A macOS binary is provided in the [latest
release](https://github.com/pnnl/fast.ssgsea/releases). Users looking to
build and install the development version of `fast.ssgsea` must have the
Xcode developer tools from Apple. See <https://mac.r-project.org/tools/>
for instructions.

### Windows

No Windows binary is available, so
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) must be
installed to compile C and C++ code. Then, the development version of
`fast.ssgsea` can be installed with the code below.

### Linux

Most Linux distributions come pre-packaged with tools to compile C and
C++ code, so no extra work is needed. Users can install the development
version of `fast.ssgsea` on Linux by running the code below.

### Install

The development version of `fast.ssgsea` can be installed with either of
the following

``` r
# install.packages("pak")
pak::pak("pnnl/fast.ssgsea")
```

``` r
# install.packages("renv")
renv::install("pnnl/fast.ssgsea")
```

## Usage

### Simulate Data

We will simulate a vector of 10,000 signed gene-level statistics. We
will also simulate 20,000 gene sets by randomly sampling between 5 and
1,000 genes.

``` r
n_genes <- 1e4L # number of genes
genes <- paste0("gene", seq_len(n_genes))

# Simulate named vector of gene-level values
set.seed(9001L)
stats <- rnorm(n = n_genes)
names(stats) <- genes

# Simulate list of gene sets
n_sets <- 2e4L
min_size <- 5L
max_size <- 500L
set_sizes <- rep(max_size:min_size, length.out = n_sets)

gene_sets <- lapply(seq_len(n_sets), function(i) {
  set.seed(i)
  sample(x = genes, size = set_sizes[i])
})
names(gene_sets) <- paste0("set", seq_along(gene_sets))
```

### Runtime and Results

This shows the runtime of `fast_ssgsea` on an AMD Ryzen 5 7600X CPU with
a clock speed of 4.7 GHz. A total of 1 million permutations were used to
calculate P-values and normalized enrichment scores (NES).

``` r
library(fast.ssgsea)

# Runtime (in seconds)
system.time({
  res <- fast_ssgsea(
    stats = stats,
    gene_sets = gene_sets,
    alpha = 1,
    nperm = 1e6L,
    min_size = min_size,
    max_size = max_size,
    seed = 0L
  )
})
```

    ##    user  system elapsed 
    ##   4.181   0.040   4.146

``` r
str(res)
```

    ## 'data.frame':    20000 obs. of  8 variables:
    ##  $ set         : chr  "set15224" "set9014" "set14650" "set7155" ...
    ##  $ set_size    : int  157 415 235 290 62 439 455 389 280 27 ...
    ##  $ ES          : num  -1681 -952 -1288 -1116 -2272 ...
    ##  $ NES         : num  -5.11 -4.74 -4.8 -4.63 -4.31 ...
    ##  $ n_same_sign : int  486561 478962 483258 482342 491895 522137 477630 520005 517360 494438 ...
    ##  $ n_as_extreme: int  24 62 66 124 149 159 147 194 199 205 ...
    ##  $ p_value     : num  5.14e-05 1.32e-04 1.39e-04 2.59e-04 3.05e-04 ...
    ##  $ adj_p_value : num  0.74 0.74 0.74 0.74 0.74 ...

### Session Information

``` r
print(sessionInfo(), locale = FALSE, tzone = FALSE)
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Linux Mint 22.3
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dqrng_0.4.1            fast.ssgsea_0.1.0.9035
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     collapse_2.1.7    fastmap_1.2.0     xfun_0.57        
    ##  [5] parallel_4.6.1    knitr_1.51        htmltools_0.5.9   rmarkdown_2.31   
    ##  [9] cli_3.6.6         data.table_1.18.4 compiler_4.6.1    rstudioapi_0.18.0
    ## [13] tools_4.6.1       evaluate_1.0.5    Rcpp_1.1.1-1.1    yaml_2.3.12      
    ## [17] otel_0.2.0        rlang_1.2.0

## Benchmarking

Benchmarking was performed on a desktop computer with an AMD Ryzen 5
7600X CPU (4.7 GHz), single threaded, to measure the runtime of
fast-ssGSEA (`fast.ssgsea::fast_ssgsea`) and FGSEA-simple
(`fgsea::fgseaSimple`). Different combinations of the number of gene
sets, maximum gene set size, and the number of permutations ($\pi$) were
tested in a random order (3 replicates each) to minimize the influence
of previous runs. The R scripts and data are available in the
simulation/ directory.

### fast-ssGSEA

<div class="figure" style="text-align: center">

<img src="./man/figures/README-figure-1.png" alt="Runtime of fast_ssgsea with 10,000, 100,000, or 1,000,000 permutations." width="720" />
<p class="caption">

Runtime of fast_ssgsea with 10,000, 100,000, or 1,000,000 permutations.
</p>

</div>

### FGSEA-simple

Like fast-ssGSEA, FGSEA-simple relies purely on the number of
permutations to calculate p-values, which limits how small they can
become. While FGSEA-simple is meant to be run with a smaller number of
permutations and followed up by FGSEA-multilevel (the method capable of
calculating arbitrarily small p-values) ([Korotkevich et al.
2021](#ref-korotkevich-fast-2021)), these results serve to illustrate
the extreme difference in runtime between the two approaches. This
difference is largely the result of changes to how the ES is defined.

<div class="figure" style="text-align: center">

<img src="./man/figures/README-figure-2.png" alt="Runtime of fgsea::fgseaSimple with 10,000, 100,000, or 1,000,000 permutations." width="720" />
<p class="caption">

Runtime of fgsea::fgseaSimple with 10,000, 100,000, or 1,000,000
permutations.
</p>

</div>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-barbie-systematic-2009" class="csl-entry">

Barbie, David A., Pablo Tamayo, Jesse S. Boehm, et al. 2009. “Systematic
RNA Interference Reveals That Oncogenic KRAS-Driven Cancers Require
TBK1.” *Nature* 462 (7269): 108–12.
<https://doi.org/10.1038/nature08460>.

</div>

<div id="ref-korotkevich-fast-2021" class="csl-entry">

Korotkevich, Gennady, Vladimir Sukhov, Nikolay Budin, Boris Shpak, Maxim
N. Artyomov, and Alexey Sergushichev. 2021. *Fast Gene Set Enrichment
Analysis*. bioRxiv. <https://doi.org/10.1101/060012>.

</div>

<div id="ref-krug-curated-2019" class="csl-entry">

Krug, Karsten, Philipp Mertins, Bin Zhang, et al. 2019. “A Curated
Resource for Phosphosite-Specific Signature Analysis.” *Molecular &
Cellular Proteomics* 18 (3): 576–93.
<https://doi.org/10.1074/mcp.TIR118.000943>.

</div>

<div id="ref-R-core-team" class="csl-entry">

R Core Team. 2026. *R: A Language and Environment for Statistical
Computing*. R Foundation for Statistical Computing.
<https://doi.org/10.32614/R.manuals>.

</div>

<div id="ref-subramanian-gene-2005" class="csl-entry">

Subramanian, Aravind, Pablo Tamayo, Vamsi K. Mootha, et al. 2005. “Gene
Set Enrichment Analysis: A Knowledge-Based Approach for Interpreting
Genome-Wide Expression Profiles.” *Proceedings of the National Academy
of Sciences* 102 (43): 15545–50.
<https://doi.org/10.1073/pnas.0506580102>.

</div>

</div>
