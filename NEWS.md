# fast.ssgsea (development version)

- Fixed Windows installation error ([#1](https://github.com/pnnl/fast.ssgsea/issues/1)).
- Greatly reduced runtime of permutation tests.
- Added `read_gmt` function, which reads a list of gene sets from a Gene Matrix Transposed (GMT) file.
- Added `alternative` parameter to `fast_ssgsea` to perform one-sided hypothesis tests.
- Added `max_size` parameter to `fast_ssgsea` to limit the maximum size of sets that will be tested.
- Updated runtime data and figures in simulation/.


# fast.ssgsea 0.1.0

- Initial release (build fails on Windows).
