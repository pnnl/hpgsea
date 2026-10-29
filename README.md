
- [hpgsea](#hpgsea)
  - [Overview](#overview)
  - [Installation](#installation)
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

The **hpgsea** R package provides a high-performance variant of
gene-permutation gene set enrichment analysis (GSEA) ([Subramanian et
al. 2005](#ref-subramanian-gene-2005)) that serves as a replacement for
the Broad Institute’s
[ssGSEA2.0](https://github.com/broadinstitute/ssGSEA2.0) repository (R
package version by Nicole Gay: <https://github.com/nicolerg/ssGSEA2>).

Benchmarks showed that **hpgsea** was at least 15.21-fold faster than
ssGSEA2.0. Furthermore, its speed, statistical power, and false positive
rates were comparable to those of the popular
[**fgsea**](https://bioconductor.org/packages//release/bioc/html/fgsea.html)
R package ([Korotkevich et al. 2021](#ref-korotkevich-fast-2021)) across
thousands of simulated datasets and 605 mouse and human datasets.
Although **hpgsea** can not compute arbitrarily small P-values like
**fgsea**, it has the advantage of being able to test bidirectional
databases—those where each biomolecule in a set has an expected
direction of change from prior experiments—like the Post-Translational
Modification Signatures Database (PTMsigDB). (Note that the combination
of **hpgsea** with PTMsigDB is known as PTM Signature Enrichment
Analysis, or PTM-SEA ([Krug et al. 2019](#ref-krug-curated-2019)).)

Moreover, while ssGSEA2.0’s false positive rate is inflated because it
unintentionally permutes the values of the genes, rather than their
labels (see [this
issue](https://github.com/broadinstitute/ssGSEA2.0/issues/45) for
details), this is not the case with **hpgsea**.

All code to reproduce the results are provided in the manuscript/
directory. The data and results will be available for download from
Zenodo in the near future.

## Installation

**macOS:** macOS users must install the [Xcode developer
tools](https://mac.r-project.org/tools/) from Apple to compile C and C++
code.

**Windows:** Windows users must install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) to compile C
and C++ code.

**Linux:** No extra work should be necessary to install **hpgsea** on
Linux, as most distributions already include a C/C++ compiler.

The development version of **hpgsea**, which requires R \>= 4.0.0, can
then be installed with

``` r
# install.packages("remotes")
remotes::install_github("pnnl/hpgsea")
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

# Simulate named list of gene sets
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

This shows the runtime of `hpgsea` on an AMD Ryzen 5 7600X (Zen 4) CPU
with its clock speed fixed at 4.7 GHz. A total of 1 million permutations
were used to calculate P-values and normalized enrichment scores (NES).

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
    ##   3.416   0.041   3.376

``` r
str(res)
```

    ## 'data.frame':    20000 obs. of  8 variables:
    ##  $ set         : chr  "set15224" "set14650" "set9014" "set7155" ...
    ##  $ set_size    : int  157 235 415 290 439 455 389 62 27 365 ...
    ##  $ ES          : num  -1681 -1288 -952 -1116 909 ...
    ##  $ NES         : num  -5.1 -4.8 -4.74 -4.64 4.43 ...
    ##  $ n_same_sign : int  486578 483605 478186 482201 522373 477393 520849 491201 493988 520384 ...
    ##  $ n_as_extreme: int  12 52 69 103 152 144 160 167 190 202 ...
    ##  $ p_value     : num  2.67e-05 1.10e-04 1.46e-04 2.16e-04 2.93e-04 ...
    ##  $ adj_p_value : num  0.534 0.759 0.759 0.759 0.759 ...

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
    ##  [7] parallel_4.6.1      knitr_1.52          htmltools_0.5.9    
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

<div id="ref-subramanian-gene-2005" class="csl-entry">

Subramanian, Aravind, Pablo Tamayo, Vamsi K. Mootha, et al. 2005. “Gene
Set Enrichment Analysis: A Knowledge-Based Approach for Interpreting
Genome-Wide Expression Profiles.” *Proceedings of the National Academy
of Sciences* 102 (43): 15545–50.
<https://doi.org/10.1073/pnas.0506580102>.

</div>

</div>
