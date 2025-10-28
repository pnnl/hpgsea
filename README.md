
- [fast.ssgsea](#fastssgsea)
  - [Installation](#installation)
    - [macOS](#macos)
    - [Windows](#windows)
    - [Install](#install)
  - [Usage](#usage)
    - [Simulate Data](#simulate-data)
    - [Runtime and Results](#runtime-and-results)
    - [Session Information](#session-information)
  - [Performance](#performance)
  - [References](#references)

# fast.ssgsea

<!-- badges: start -->

[![R-CMD-check](https://github.com/pnnl/fast.ssgsea/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pnnl/fast.ssgsea/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/394311897.svg)](https://doi.org/10.5281/zenodo.16783102)
<!-- badges: end -->

`fast.ssgsea` is an R package ([R Core Team 2024](#ref-R-core-team)) for
fast Single-Sample Gene Set Enrichment Analysis (ssGSEA) and
Post-Translational Modification Signature Enrichment Analysis (PTM-SEA)
([Barbie et al. 2009](#ref-barbie-systematic-2009); [Krug et al.
2019](#ref-krug-curated-2019)).

The primary function, `fast_ssgsea`, accepts a numeric matrix with genes
or other molecules as rows and either samples, contrasts, or some other
meaningful representation of the data as columns. A named list of gene
sets (more generally, molecular signatures) is also required. Other
arguments control the behavior of ssGSEA/PTM-SEA, and they are described
in the function documentation.

The package also contains a `read_gmt` function, which reads a Gene
Matrix Transposed (GMT) file to construct a named list of gene sets for
use with `fast_ssgsea`.

## Installation

R version 4.0.0 or greater is required to install `fast.ssgsea`.

It may be possible to get `fast.ssgsea` to work with older versions of R
by cloning the repository, changing the minimum R version in the
DESCRIPTION (e.g., to `>= 3.6.0`), and rebuilding the package, but users
should do so at their own risk.

### macOS

A macOS binary is provided in the [latest
release](https://github.com/pnnl/fast.ssgsea/releases). Users looking to
build and install the development version of `fast.ssgsea` must have the
Xcode developer tools from Apple and a FORTRAN compiler installed. See
<https://mac.r-project.org/tools/> for instructions.

### Windows

No Windows binary is available, so
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) must be
installed to compile C++ code.

### Install

The development version of `fast.ssgsea` can be installed with

``` r
# install.packages("pak")
pak::pak("pnnl/fast.ssgsea")
```

## Usage

### Simulate Data

We will simulate a matrix with 10,000 genes as rows and 100 samples as
columns. Then, we generate 20,000 gene sets by randomly sampling between
10 and 500 genes from the matrix row names.

``` r
n_genes <- 10000L # number of genes
n_samples <- 100L # number of samples
genes <- paste0("gene", seq_len(n_genes))
samples <- paste0("sample", seq_len(n_samples))

## Simulate matrix of sample gene expression values
set.seed(9001L)
X <- matrix(
  data = rnorm(n = n_genes * n_samples),
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(genes, samples)
)

## Simulate list of gene sets
n_sets <- 20000L # number of gene sets
min_size <- 5L # size of smallest gene set
max_size <- 1000L # size of largest gene set

size_range <- max_size - min_size + 1L
n_reps <- ceiling(n_sets / size_range)
set_sizes <- rep(max_size:min_size, times = n_reps)[seq_len(n_sets)]

gene_sets <- lapply(seq_len(n_sets), function(i) {
  set.seed(i)
  sample(x = genes, size = set_sizes[i])
})
names(gene_sets) <- paste0("set", seq_along(gene_sets))
```

### Runtime and Results

This shows the runtime of `fast_ssgsea` running on an AMD Ryzen 5 7600X
CPU with a clock speed of 4.7 GHz.

``` r
library(fast.ssgsea)

# Runtime (elapsed time)
system.time({
  res <- fast_ssgsea(
    X = X,
    gene_sets = gene_sets,
    alpha = 1,
    nperm = 1000L,
    batch_size = 1000L,
    adjust_globally = FALSE,
    min_size = min_size,
    sort = TRUE,
    seed = 0L
  )
})
```

    ##    user  system elapsed 
    ##  15.572   1.352   9.120

``` r
str(res)
```

    ## 'data.frame':    2000000 obs. of  9 variables:
    ##  $ sample      : Factor w/ 100 levels "sample1","sample2",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ set         : chr  "set4576" "set12526" "set11427" "set9645" ...
    ##  $ set_size    : int  409 427 530 320 320 977 519 517 511 841 ...
    ##  $ ES          : num  929 861 693 1043 898 ...
    ##  $ NES         : num  4.4 4.13 3.72 4.22 3.64 ...
    ##  $ n_same_sign : int  544 539 536 534 534 525 521 521 521 520 ...
    ##  $ n_as_extreme: int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ p_value     : num  0.00183 0.00185 0.00186 0.00187 0.00187 ...
    ##  $ adj_p_value : num  0.838 0.838 0.838 0.838 0.838 ...

### Session Information

``` r
print(sessionInfo(), locale = FALSE, tzone = FALSE)
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Linux Mint 22.1
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fast.ssgsea_0.1.0.9017
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] dqrng_0.4.1            digest_0.6.37          RcppArmadillo_15.0.2-2
    ##  [4] fastmap_1.2.0          xfun_0.53              Matrix_1.7-4          
    ##  [7] lattice_0.22-7         knitr_1.50             htmltools_0.5.8.1     
    ## [10] rmarkdown_2.29         cli_3.6.5              grid_4.5.1            
    ## [13] data.table_1.17.8      compiler_4.5.1         rstudioapi_0.17.1     
    ## [16] tools_4.5.1            evaluate_1.0.5         Rcpp_1.1.0            
    ## [19] yaml_2.3.10            rlang_1.1.6

## Performance

The `fast.ssgsea` R package utilizes linear algebra and ideas from Fast
Gene Set Enrichment Analysis ([Korotkevich et al.
2021](#ref-korotkevich-fast-2021)) to greatly reduce the runtime of gene
permutation GSEA and PTM-SEA.

Tests were performed on a desktop computer with an AMD Ryzen 5 7600X CPU
(6 cores, 12 threads) at 4.7 GHz. Different combinations of the number
of samples, gene sets, maximum gene set size, number of permutations,
and value of the $\alpha$ parameter (the weighting exponent) were tested
in a random order (3 replicates each) to minimize the influence of
previous runs.

<div class="figure" style="text-align: center">

<img src="./man/figures/README-figure-1.png" alt="Runtime of fast_ssgsea with A) 1,000 or B) 10,000 permutations." width="749" />
<p class="caption">

Runtime of fast_ssgsea with A) 1,000 or B) 10,000 permutations.
</p>

</div>

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-barbie-systematic-2009" class="csl-entry">

Barbie, David A., Pablo Tamayo, Jesse S. Boehm, So Young Kim, Susan E.
Moody, Ian F. Dunn, Anna C. Schinzel, et al. 2009. “Systematic RNA
Interference Reveals That Oncogenic KRAS-Driven Cancers Require TBK1.”
*Nature* 462 (7269): 108–12. <https://doi.org/10.1038/nature08460>.

</div>

<div id="ref-korotkevich-fast-2021" class="csl-entry">

Korotkevich, Gennady, Vladimir Sukhov, Nikolay Budin, Boris Shpak, Maxim
N. Artyomov, and Alexey Sergushichev. 2021. “Fast Gene Set Enrichment
Analysis.” bioRxiv. <https://doi.org/10.1101/060012>.

</div>

<div id="ref-krug-curated-2019" class="csl-entry">

Krug, Karsten, Philipp Mertins, Bin Zhang, Peter Hornbeck, Rajesh Raju,
Rushdy Ahmad, Matthew Szucs, et al. 2019. “A Curated Resource for
Phosphosite-Specific Signature Analysis.” *Molecular & Cellular
Proteomics* 18 (3): 576–93. <https://doi.org/10.1074/mcp.TIR118.000943>.

</div>

<div id="ref-R-core-team" class="csl-entry">

R Core Team. 2024. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

</div>
