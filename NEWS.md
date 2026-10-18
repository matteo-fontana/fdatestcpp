# fdatestcpp 0.1.0

* First release.
* Interval-wise (IWT), threshold-wise (TWT) and global p-value adjustment on
  precomputed observed and resampled pointwise statistics: `adjust_iwt()`,
  `adjust_twt()`, `adjust_global()`, `pvalues_pointwise()`.
* One- and two-sample resampling front ends `perm_stats1()` and
  `perm_stats2()`: sign flips, paired swaps and label permutations, optionally
  constant within blocks of observations, or user-supplied resamples.
* `iwt1()`, `iwt2()`, `twt2()` and `global2()` with the interface of the
  functions of the same names in the development version of 'fdatest'
  (<https://github.com/permaverse/fdatest>; the CRAN release 2.1.1 of
  'fdatest' does not have them), plus `twt1()` and `global1()`, for curves
  given as matrices or as 'fda' functional data objects.
* `plot()` and `autoplot()` methods for the results of these functions:
  'ggplot2' graphics after the IWTgg plotting scheme, with the curves and the
  adjusted and unadjusted p-value functions in two panels and the part of the
  domain selected by the test shaded (on grids of at most 50 points the
  p-value functions are drawn as steps over the shaded cells). The results
  have classes `c("fdatestcpp_fos", "fos")` and `c("fdatestcpp_fts", "fts")`.
