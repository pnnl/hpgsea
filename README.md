
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
<!-- badges: end -->

`fast.ssgsea` is an R package ([R Core Team 2024](#ref-R-core-team)) for
High-Performance Gene Set Enrichment Analysis (HP-GSEA). It is also
capable of performing Post-Translational Modification Signature
Enrichment Analysis (PTM-SEA) ([Subramanian et al.
2005](#ref-subramanian-gene-2005); [Krug et al.
2019](#ref-krug-curated-2019)).

The primary function, `fast_ssgsea`, accepts a numeric matrix with genes
or other molecules as rows and either samples, contrasts, or some other
meaningful representation of the data as columns. A named list of gene
sets (more generally, molecular signatures) is also required. Other
arguments control the behavior of HP-GSEA/PTM-SEA, and they are
described in the function documentation.

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

We will simulate a matrix with 10,000 genes as rows and one column.
Then, we generate 20,000 gene sets by randomly sampling between 5 and
1,000 genes.

``` r
n_genes <- 10000L # number of genes
n_samples <- 1L # number of samples (>= 1)
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

This shows the runtime of `fast_ssgsea` on an AMD Ryzen 5 7600X CPU with
a clock speed of 4.7 GHz. A total of 100,000 permutations were used to
calculate P-values and normalized enrichment scores (NES).

``` r
library(fast.ssgsea)

# Runtime (elapsed time)
system.time({
  res <- fast_ssgsea(
    X = X,
    gene_sets = gene_sets,
    alpha = 1,
    nperm = 1e5L, # default is 1000
    min_size = min_size,
    seed = 0L
  )
})
```

    ##    user  system elapsed 
    ##   5.919   0.906   6.403

``` r
str(res)
```

    ## 'data.frame':    20000 obs. of  9 variables:
    ##  $ sample      : Factor w/ 1 level "sample1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ set         : chr  "set18791" "set16136" "set19084" "set2830" ...
    ##  $ set_size    : int  138 801 841 163 706 749 450 87 161 761 ...
    ##  $ ES          : num  -1866 709 698 1584 759 ...
    ##  $ NES         : num  -5.3 4.65 4.68 4.76 4.68 ...
    ##  $ n_same_sign : int  49042 52788 52782 50951 52785 47193 47813 50722 48979 47243 ...
    ##  $ n_as_extreme: int  1 8 8 9 11 10 13 14 18 20 ...
    ##  $ p_value     : num  4.08e-05 1.70e-04 1.71e-04 1.96e-04 2.27e-04 ...
    ##  $ adj_p_value : num  0.739 0.739 0.739 0.739 0.739 ...

``` r
head(res, 10L)
```

    ##     sample      set set_size         ES       NES n_same_sign n_as_extreme
    ## 1  sample1 set18791      138 -1865.9539 -5.301993       49042            1
    ## 2  sample1 set16136      801   709.4930  4.647777       52788            8
    ## 3  sample1 set19084      841   697.9020  4.677156       52782            8
    ## 4  sample1  set2830      163  1584.3635  4.761980       50951            9
    ## 5  sample1 set18223      706   759.1365  4.680156       52785           11
    ## 6  sample1 set17184      749  -668.8030 -4.503202       47193           10
    ## 7  sample1  set8519      450  -859.1930 -4.454985       47813           13
    ## 8  sample1 set10874       87  1927.4191  4.267786       50722           14
    ## 9  sample1  set2832      161 -1380.4446 -4.245356       48979           18
    ## 10 sample1  set3228      761  -633.4807 -4.310214       47243           20
    ##         p_value adj_p_value
    ## 1  4.078054e-05   0.7393096
    ## 2  1.704901e-04   0.7393096
    ## 3  1.705094e-04   0.7393096
    ## 4  1.962631e-04   0.7393096
    ## 5  2.273330e-04   0.7393096
    ## 6  2.330805e-04   0.7393096
    ## 7  2.928013e-04   0.7393096
    ## 8  2.957238e-04   0.7393096
    ## 9  3.879134e-04   0.7984999
    ## 10 4.445009e-04   0.7984999

### Session Information

``` r
print(sessionInfo(), locale = FALSE, tzone = FALSE)
```

    ## R version 4.5.2 (2025-10-31)
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
    ## [1] fast.ssgsea_0.1.0.9018
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] dqrng_0.4.1            digest_0.6.37          RcppArmadillo_15.0.2-2
    ##  [4] fastmap_1.2.0          xfun_0.53              Matrix_1.7-4          
    ##  [7] lattice_0.22-7         knitr_1.50             htmltools_0.5.8.1     
    ## [10] rmarkdown_2.29         cli_3.6.5              grid_4.5.2            
    ## [13] data.table_1.17.8      compiler_4.5.2         rstudioapi_0.17.1     
    ## [16] tools_4.5.2            evaluate_1.0.5         Rcpp_1.1.0            
    ## [19] yaml_2.3.10            rlang_1.1.6

## Performance

The `fast.ssgsea` R package utilizes linear algebra and ideas from Fast
Gene Set Enrichment Analysis ([Korotkevich et al.
2021](#ref-korotkevich-fast-2021)) to greatly reduce the runtime.

Tests were performed on a desktop computer with an AMD Ryzen 5 7600X CPU
(6 cores, 12 threads) at 4.7 GHz. Different combinations of the number
of gene sets, maximum gene set size, number of permutations, and value
of the $\alpha$ parameter (the weighting exponent) were tested in a
random order (3 replicates each) to minimize the influence of previous
runs.

<div class="figure" style="text-align: center">

<img src="./man/figures/README-figure-1.png" alt="Runtime of fast_ssgsea with A) 10,000, B) 100,000, or C) 1,000,000 permutations." width="648" />
<p class="caption">

Runtime of fast_ssgsea with A) 10,000, B) 100,000, or C) 1,000,000
permutations.
</p>

</div>

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

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

<div id="ref-subramanian-gene-2005" class="csl-entry">

Subramanian, Aravind, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee,
Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, et al. 2005.
“Gene Set Enrichment Analysis: A Knowledge-Based Approach for
Interpreting Genome-Wide Expression Profiles.” *Proceedings of the
National Academy of Sciences* 102 (43): 15545–50.
<https://doi.org/10.1073/pnas.0506580102>.

</div>

</div>
