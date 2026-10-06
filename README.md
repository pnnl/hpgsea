
- [hpgsea](#hpgsea)
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
  - [References](#references)

# hpgsea

<!-- badges: start -->

[![R-CMD-check](https://github.com/pnnl/hpgsea/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pnnl/hpgsea/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`hpgsea` is an R package ([R Core Team 2026](#ref-R-core-team)) for a
highly optimized variant of pre-ranked Gene Set Enrichment Analysis
(GSEA) ([Subramanian et al. 2005](#ref-subramanian-gene-2005)). Unlike
standard GSEA, HPGSEA is capable of testing gene sets where each gene
has an expected direction of change (up- or down-regulation; indicated
by appending a “;u” or “;d” to the end of every gene in a set) from a
prior experiment.

HPGSEA is based on Post-Translational Modification Signature Enrichment
Analysis (PTM-SEA) ([Krug et al. 2019](#ref-krug-curated-2019)), and it
borrows optimization techniques from the simple implementation of Fast
Gene Set Enrichment Analysis (FGSEA-simple) ([Korotkevich et al.
2021](#ref-korotkevich-fast-2021)).

The primary function, `hpgsea`, accepts a vector of signed statistics
with genes or other molecules as names. The values must be approximately
symmetric around zero, with more extreme values indicating greater
importance. A named list of gene sets (more generally, molecular
signatures) is also required. Other arguments control the behavior of
HPGSEA, and they are described in the documentation.

The package also contains a `read_gmt` function, which reads a Gene
Matrix Transposed (GMT) file to construct a named list of gene sets for
use with `hpgsea`.

## Installation

R version 4.0.0 or greater is required to install `hpgsea`.

### macOS

A macOS binary is provided in the [latest
release](https://github.com/pnnl/hpgsea/releases). Users looking to
build and install the development version of `hpgsea` must have the
Xcode developer tools from Apple. See <https://mac.r-project.org/tools/>
for instructions.

### Windows

No Windows binary is available, so
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) must be
installed to compile C and C++ code. Then, the development version of
`hpgsea` can be installed with the code below.

### Linux

Most Linux distributions come pre-packaged with tools to compile C and
C++ code, so no extra work is needed. Users can install the development
version of `hpgsea` on Linux by running the code below.

### Install

The development version of `hpgsea` can be installed with either of the
following

``` r
# install.packages("pak")
pak::pak("pnnl/hpgsea")
```

``` r
# install.packages("renv")
renv::install("pnnl/hpgsea")
```

## Usage

### Simulate Data

We will simulate a vector of 10,000 signed gene-level statistics and a
list of 20,000 gene sets by randomly sampling between 5 and 500 genes.

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

This shows the runtime of `hpgsea` on an AMD Ryzen 5 7600X CPU with a
clock speed of 4.7 GHz. A total of 1 million permutations were used to
calculate P-values and normalized enrichment scores (NES).

``` r
library(hpgsea)

# Runtime (in seconds)
set.seed(0L)
system.time({
  res <- hpgsea(
    stats = stats,
    gene_sets = gene_sets,
    alpha = 1,
    nperm = 1e6L,
    min_size = min_size,
    max_size = max_size
  )
})
```

    ##    user  system elapsed 
    ##   3.499   0.039   3.457

``` r
str(res)
```

    ## 'data.frame':    20000 obs. of  8 variables:
    ##  $ set         : chr  "set15224" "set14650" "set9014" "set7155" ...
    ##  $ set_size    : int  157 235 415 290 439 62 455 389 280 365 ...
    ##  $ ES          : num  -1681 -1288 -952 -1116 909 ...
    ##  $ NES         : num  -5.1 -4.8 -4.74 -4.63 4.43 ...
    ##  $ n_same_sign : int  486027 483686 478465 482064 522291 490985 477346 520334 517260 520177 ...
    ##  $ n_as_extreme: int  18 65 81 129 141 145 154 175 179 201 ...
    ##  $ p_value     : num  3.91e-05 1.36e-04 1.71e-04 2.70e-04 2.72e-04 ...
    ##  $ adj_p_value : num  0.698 0.698 0.698 0.698 0.698 ...

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
    ## [1] hpgsea_0.1.0.9037
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39       dqrng_0.4.1         collapse_2.1.8     
    ##  [4] fastmap_1.2.0       xfun_0.60           BH_1.90.0-1        
    ##  [7] parallel_4.6.1      knitr_1.51          htmltools_0.5.9    
    ## [10] rmarkdown_2.31      cli_3.6.6           data.table_1.18.6.1
    ## [13] compiler_4.6.1      rstudioapi_0.19.0   tools_4.6.1        
    ## [16] evaluate_1.0.5      Rcpp_1.1.2          yaml_2.3.12        
    ## [19] otel_0.2.0          rlang_1.3.0

## References

<div id="refs" class="references csl-bib-body hanging-indent">

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
