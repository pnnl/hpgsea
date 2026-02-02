
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
    - [fast-ssGSEA](#fast-ssgsea)
    - [FGSEA-simple](#fgsea-simple)
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
    nperm = 1e5L,
    min_size = min_size,
    seed = 0L
  )
})
```

    ##    user  system elapsed 
    ##   2.245   0.377   2.376

``` r
str(res)
```

    ## 'data.frame':    20000 obs. of  9 variables:
    ##  $ sample      : Factor w/ 1 level "sample1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ set         : chr  "set18791" "set2830" "set19084" "set18223" ...
    ##  $ set_size    : int  138 163 841 706 801 87 503 409 320 450 ...
    ##  $ ES          : num  -1866 1584 698 759 709 ...
    ##  $ NES         : num  -5.34 4.78 4.66 4.67 4.62 ...
    ##  $ n_same_sign : int  49235 51107 52907 52784 52813 50461 52351 51847 51728 47859 ...
    ##  $ n_as_extreme: int  1 3 8 9 12 12 16 19 19 19 ...
    ##  $ p_value     : num  4.06e-05 7.83e-05 1.70e-04 1.89e-04 2.46e-04 ...
    ##  $ adj_p_value : num  0.783 0.783 0.836 0.836 0.836 ...

### Session Information

``` r
print(sessionInfo(), locale = FALSE, tzone = FALSE)
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Linux Mint 22.1
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dqrng_0.4.1            fast.ssgsea_0.1.0.9025
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.37          RcppArmadillo_15.0.2-2 collapse_2.1.3        
    ##  [4] fastmap_1.2.0          xfun_0.53              Matrix_1.7-4          
    ##  [7] lattice_0.22-7         parallel_4.5.2         knitr_1.50            
    ## [10] htmltools_0.5.8.1      rmarkdown_2.29         cli_3.6.5             
    ## [13] grid_4.5.2             data.table_1.17.8      compiler_4.5.2        
    ## [16] rstudioapi_0.17.1      tools_4.5.2            evaluate_1.0.5        
    ## [19] Rcpp_1.1.0             yaml_2.3.10            rlang_1.1.6

## Performance

### fast-ssGSEA

The `fast.ssgsea` R package utilizes linear algebra and ideas from Fast
Gene Set Enrichment Analysis (FGSEA) ([Korotkevich et al.
2021](#ref-korotkevich-fast-2021)) to greatly reduce the runtime.

Tests were performed on a desktop computer with an AMD Ryzen 5 7600X CPU
running at 4.7 GHz, single threaded. Different combinations of the
number of gene sets, maximum gene set size, and the number of
permutations were tested in a random order (3 replicates each) to
minimize the influence of previous runs. The R scripts and data are
available in the simulation/ folder.

<div class="figure" style="text-align: center">

<img src="./man/figures/README-figure-1.png" alt="Runtime of fast_ssgsea with A) 10,000, B) 100,000, or C) 1,000,000 permutations." width="864" />
<p class="caption">

Runtime of fast_ssgsea with A) 10,000, B) 100,000, or C) 1,000,000
permutations.
</p>

</div>

### FGSEA-simple

The same tests were also carried out using the simple implementation of
FGSEA (`fgsea::fgseaSimple`). Like fast-ssGSEA, FGSEA-simple relies
purely on the number of permutations to calculate p-values, which limits
how small they can become. While FGSEA-simple is instead meant to be run
with a smaller number of permutations and followed up with
FGSEA-multilevel (the method capable of calculating arbitrarily small
p-values), these results serve to illustrate the extreme difference in
runtime between the two approaches. This difference is largely the
result of changes to how the ES is defined.

<div class="figure" style="text-align: center">

<img src="./man/figures/README-figure-2.png" alt="Runtime of fgsea::fgseaSimple with A) 10,000, B) 100,000, or C) 1,000,000 permutations." width="864" />
<p class="caption">

Runtime of fgsea::fgseaSimple with A) 10,000, B) 100,000, or C)
1,000,000 permutations.
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
