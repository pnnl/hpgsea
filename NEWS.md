# fast.ssgsea (development version)

## BREAKING CHANGES
- Renamed the `X` parameter of `fast_ssgsea` to `stats`. Parameter `stats` accepts a named numeric vector, rather than a numeric matrix with dimension names. This also led to the removal of the `adjust.globally` parameter.
- Removed the `batch.size` parameter from `fast_ssgsea`.

## ENHANCEMENTS
- Greatly reduced runtime of permutation tests, especially when testing directional gene sets.
- Increased default `nperm` to 100,000 in `fast_ssgsea` and increased the permutation limit from 1 million to 2 billion, which is close to `.Machine$integer.max`.
- Greatly reduced memory usage and slightly improved runtime by removing the need for incidence matrices. This also led to the removal of the _Matrix_ and _RcppArmadillo_ packages from Imports.
- Added parameter `alternative` to `fast_ssgsea` to perform one-sided hypothesis tests.
- Added parameter `max_size` to `fast_ssgsea` to limit the maximum size of sets that will be tested.
- Added function `read_gmt`, which reads a named list of gene sets from a Gene Matrix Transposed (GMT) file.

## BUGFIXES
- Fixed Windows installation error ([#1](https://github.com/pnnl/fast.ssgsea/issues/1)).
- Directional gene sets are now allowed to consist entirely of up-regulated or down-regulated genes.
- Permutation enrichment scores for down-regulated genes are now calculated so they avoid overlap with the genes selected for the up-regulated permutation enrichment scores.

## MISC
- Added the _collapse_ package to Imports.
- Updated runtime data and figures in simulation/.


# fast.ssgsea 0.1.0

- Initial release (build fails on Windows).
