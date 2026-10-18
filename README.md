# fdatestcpp

Fast interval-wise and threshold-wise tests for functional data in R.

fdatestcpp implements in C++ the one- and two-sample permutation tests of the
[fdatest](https://github.com/permaverse/fdatest) package by Alessia Pini,
Aymeric Stamm and Simone Vantini:

- the **interval-wise testing procedure** (IWT; Pini and Vantini, 2017),
- the **threshold-wise testing procedure** (TWT; Abramowicz et al., 2023),
- the **global test** and unadjusted pointwise p-values.

These tests find the parts of the domain where two mean functions differ (or
where a mean differs from a reference), with control of the family-wise error
rate on intervals (IWT) or on sub- and super-level sets (TWT).

- **Same interface.** `iwt1()`, `iwt2()`, `twt2()` and `global2()` have the
  argument names, defaults and output fields of the functions of the same
  names in fdatest, and take curves as [fda](https://CRAN.R-project.org/package=fda)
  `fd` objects or as matrices.
- **Same results.** With the same inputs and random seed, the p-values,
  adjusted p-values and IWT p-value matrices are `identical()` to fdatest's
  (see [Exactness](#exactness) for the one deliberate exception and the
  platform condition).
- **Faster.** The IWT adjustment costs O(B p^2) operations instead of
  O(B p^3) for B resamples and p grid points, and the resampling runs in C++:
  17 to 1,074 times faster than fdatest end to end in the
  [benchmarks](#speed).
- **More resampling schemes.** Sign flips and swaps constant within blocks of
  consecutive observations (for dependent curves), user-supplied sign flips
  or permutations, and an engine that adjusts p-values for any statistics you
  compute yourself.
- **Plots.** `plot()` and `autoplot()` draw the curves and the p-value
  functions with ggplot2, shading the part of the domain the test selects.

fdatestcpp follows the **development version** of fdatest (GitHub
`permaverse/fdatest`, commit `67bb41c`, April 2026). The CRAN release of
fdatest (2.1.1, 2022) has the earlier interval testing procedure (`ITP*`
functions) but not `iwt1()` and the other functions ported here; the
development version also reports version number 2.1.1.

## Installation

```r
# install.packages("remotes")
remotes::install_github("matteo-fontana/fdatestcpp")
```

The package is compiled from source: on Windows this needs
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), on macOS the Xcode
command line tools. It depends on Rcpp, fda and ggplot2.

fdatest itself is not needed. To compare results with it, install its
development version from r-universe:

```r
install.packages("fdatest",
                 repos = c("https://alessiapini.r-universe.dev", "https://cloud.r-project.org"))
```

fdatestcpp uses fdatest's function names. If both packages are attached, the
one attached last masks the other: call `fdatestcpp::iwt1()` explicitly, or
attach fdatestcpp after fdatest.

## Quick start

```r
library(fdatestcpp)
library(fda)

# 50 smooth curves on [0, 1] whose mean has a bump around t = 0.7,
# and 50 curves with mean zero
set.seed(1)
tt <- seq(0, 1, length.out = 200)
noise <- function(n) {
  phi <- cbind(sapply(1:4, function(k) cos(2 * pi * k * tt) / k),
               sapply(1:4, function(k) sin(2 * pi * k * tt) / k))
  matrix(rnorm(n * 8, sd = 0.5), n, 8) %*% t(phi)
}
y <- sweep(noise(50), 2, 0.8 * exp(-((tt - 0.7) / 0.1)^2), "+")
z <- noise(50)

# fd objects, as returned by fda::smooth.basis()
basis <- create.bspline.basis(c(0, 1), nbasis = 30)
y_fd <- smooth.basis(tt, t(y), basis)$fd
z_fd <- smooth.basis(tt, t(z), basis)$fd

# One-sample interval-wise test of H0: mean = 0 on a non-periodic domain
r1 <- iwt1(y_fd, mu = 0, n_perm = 2000, recycle = FALSE)
r1$global_pvalue                                  # 0
range(r1$abscissa[r1$adjusted_pvalues <= 0.05])   # 0.64 0.78

# Two-sample interval-wise and threshold-wise tests of equal means
r2 <- iwt2(y_fd, z_fd, n_perm = 2000, recycle = FALSE)
r3 <- twt2(y_fd, z_fd, n_perm = 2000)
range(r2$abscissa[r2$adjusted_pvalues <= 0.05])   # 0.67 0.75
range(r3$abscissa[r3$adjusted_pvalues <= 0.05])   # 0.65 0.77

# Curves and p-value functions, with the selected domain shaded
plot(r1)
plot(r3, title = "Two-sample TWT")
```

![Output of plot(r1): the curves, their mean and mu above; the adjusted and
unadjusted p-value functions below; the domain where the adjusted p-value is
below 0.05 shaded in grey, darker below 0.01.](man/figures/README-plot-iwt1.png)

![Output of plot(r3): the curves of the two samples in orange and blue with
their means; the TWT-adjusted and unadjusted p-value functions; the selected
domain shaded.](man/figures/README-plot-twt2.png)

As in fdatest, p-values are the fraction of resamples whose statistic is at
least the observed one, so they can be 0; `plus_one = TRUE` gives
`(1 + count) / (B + 1)` instead. `recycle = TRUE` (fdatest's default) treats
the domain as periodic, so that intervals wrap around its end.

More examples:

```r
# Matrices: one curve per row, one grid point per column
r <- iwt1(y, n_perm = 2000, recycle = FALSE)

# Dependent curves (e.g. a time series of curves): signs constant within
# 10 blocks of consecutive curves instead of curve by curve
r <- iwt1(y_fd, n_perm = 2000, recycle = FALSE, blocks = 10)

# Any resampling scheme: adjust p-values for statistics computed elsewhere
st <- perm_stats1(y_fd, n_perm = 2000)         # or your own t0 and t_perm
iwt <- adjust_iwt(st$t0, st$t_perm, recycle = FALSE)
twt <- adjust_twt(st$t0, st$t_perm)
```

## Functions

The package has three layers. The wrappers combine the other two.

### Wrappers

`iwt1()`, `iwt2()`, `twt2()` and `global2()` have fdatest's arguments,
defaults and output fields. `twt1()` and
`global1()` are one-sample versions not in fdatest. On top of fdatest's
arguments they take:

- `alternative` and `standardize` in one-sample tests (fdatest has them in
  two-sample tests only);
- `blocks`, `signs` and `perms` (see the resampling front end below);
- `plus_one` (see above) and `return_perm` (return the resampled statistics).

The outputs add `global_pvalue` (fdatest's `global2()` already has it), the
observed statistics `t0`, the grid `abscissa` for `fd` data and, with
`return_perm = TRUE`, the resampled statistics `t_perm`. For the IWT,
`global_pvalue` is `pvalue_matrix[1, 1]`. The outputs have fdatest's classes
`"fos"` (one sample) and `"fts"` (two samples), preceded by
`"fdatestcpp_fos"` and `"fdatestcpp_fts"` for the plot methods below.

### Plots

`plot(r)` draws a result and `autoplot(r)` returns it as a ggplot, in the
style of the [IWTgg](https://github.com/matteo-fontana/IWTgg) plotting
scheme: two panels over a common abscissa, the curves above and the
adjusted (solid) and unadjusted (dashed) p-value functions below. The part
of the domain where the adjusted p-value is below `alpha1 = 0.05` is shaded
in light grey, below `alpha2 = 0.01` in darker grey. Thick lines are the
sample means; the dashed line is the null hypothesis (`mu` in one-sample
tests, the mean of the second sample plus `mu` in two-sample tests with
`mu` not 0, explained in the caption). The two samples are drawn in
colour-blind-safe orange and blue, which
`+ ggplot2::scale_colour_manual(values = c("1" = ..., "2" = ...), name = "Sample")`
changes.

```r
plot(r, alpha1 = 0.1, alpha2 = 0.05, title = "IWT")
autoplot(r, panels = "pvalues", unadjusted = FALSE) + ggplot2::labs(caption = "B = 2000")
plot(r, xrange = c(0, 24), xlab = "Hour")   # matrix data: abscissa for the grid
ggplot2::ggsave("iwt.pdf", autoplot(r), width = 7, height = 5)
```

The abscissa is the evaluation grid for `fd` data, and `[0, 1]` for
matrices, as in fdatest; for matrices, `xrange` gives the range of the grid.
For `fd` data an `xrange` that moves the points off the evaluation grid
gives a warning (to relabel the axis, transform `r$abscissa` instead). On
coarse grids (at most 50 points) the p-value functions are drawn as steps
over the shaded cells, so that a curve crosses a level exactly where the
shading starts; `steps = TRUE` or `FALSE` forces either style. The two
panels are facets of one ggplot, so
the result can be modified and saved like any other ggplot. A complete
theme such as `theme_bw()` resets the layout; follow it with
`theme(strip.placement = "outside", legend.position = "bottom",
legend.box = "vertical", legend.location = "plot")` to keep the panel labels
on the left and the legends below.
`autoplot()` warns about arguments it does not use, such as base graphics'
`main` or `lwd`. For the global tests the adjusted p-value is the global
p-value everywhere, so the whole domain is shaded when the test rejects. fdatest's `IWTimage()`
heatmap of the p-value matrix also works on IWT results, since they keep
fdatest's classes (on one-sample results with an `fd` `mu`, set
`r$mu <- as.vector(r$mu)` first, as for fdatest's own results).

### Resampling front end

| Function | Resampling |
|---|---|
| `perm_stats1(data, mu = 0, n_perm = 1000L, dx = NULL, alternative, standardize = FALSE, blocks = NULL, signs = NULL)` | one sample: sign flips of `data - mu` |
| `perm_stats2(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE, alternative, standardize = FALSE, blocks = NULL, signs = NULL, perms = NULL)` | paired: swaps within pairs (sign flips of the differences); unpaired: label permutations |

Both return the observed pointwise statistics `t0`, the `B x p` matrix of
resampled statistics `t_perm`, the number of resampling units `n_units`, `mu`
and, for `fd` data, `abscissa`. The argument order follows fdatest's `iwt1()`
and `iwt2()` (hence `dx` before `n_perm` in `perm_stats2()`).

- `alternative`: `"two.sided"` uses the squared mean (difference);
  `"greater"` and `"less"` set it to 0 where its sign is not the
  alternative's, as fdatest does.
- `standardize = TRUE` divides by a variance recomputed on every resample:
  fdatest's pooled variance for two samples, and for one sample the variance
  of the sign-flipped curves (`colMeans(y)^2 / apply(y, 2, var)`). Where all
  curves are equal, the statistic is `Inf`, or `NaN` when its numerator is 0
  too (curves all zero there, as `fd` curves pinned to 0 at a boundary, or a
  mean of the wrong sign under a one-sided alternative). As in fdatest, a `NaN`
  statistic gives an `NA` unadjusted p-value at its point, is ignored by the
  IWT and TWT, and makes the global test's p-values `NA`.
- `blocks`: `NULL` (each curve, or pair, is a resampling unit), a number `q`
  of blocks of consecutive curves (`cut(seq_len(n), q, labels = FALSE)`), or
  a vector of block ids. Signs (swaps) are constant within a block, and each
  resample costs O(q p) instead of O(n p). With blocks, `t0` is computed from
  block sums and can differ from the unblocked `t0` in the last bits.
- `signs`: a user-supplied `B x q` matrix of -1 and 1, one column per unit
  (for a vector of block ids, column `g` is `unique(blocks)[g]`). One sample:
  -1 flips the unit's curves. Paired: -1 swaps them.
- `perms`: unpaired tests only, a `B x (n1 + n2)` matrix whose rows are
  permutations of the rows of `rbind(data1, data2)`.

### Engine

`t0` is the length-`p` vector of observed pointwise statistics and `t_perm`
the `B x p` matrix of resampled ones, with large values speaking against the
null hypothesis. They can come from `perm_stats1()`/`perm_stats2()` or from
any other resampling scheme.

| Function | Returns | fdatest counterpart |
|---|---|---|
| `pvalues_pointwise(t0, t_perm, plus_one = FALSE)` | unadjusted p-values | computed inside `*_permtest()` |
| `adjust_iwt(t0, t_perm, recycle, aggregation_strategy = c("integral", "max"), plus_one = FALSE, pvalue_matrix = TRUE)` | `unadjusted_pvalues`, `adjusted_pvalues`, `pvalue_matrix`, `global_pvalue` | `fdatest:::p_adjust_iwt()` |
| `adjust_twt(t0, t_perm, aggregation_strategy, plus_one = FALSE)` | `unadjusted_pvalues`, `adjusted_pvalues`, `global_pvalue` | `fdatest:::p_adjust_twt()` |
| `adjust_global(t0, t_perm, aggregation_strategy, plus_one = FALSE)` | `unadjusted_pvalues`, `adjusted_pvalues` (constant), `global_pvalue` | `fdatest:::p_adjust_global()` |

`adjust_iwt()` has no default for `recycle`, so that the choice between a
periodic and a non-periodic domain is explicit. Code that calls fdatest's
internal adjustment on its own statistics, such as

```r
fdatest:::p_adjust_iwt(p = p, pval = pval, t0 = t0, t_coeff = t_coeff,
                       n_perm = B, recycle = FALSE, aggregation_strategy = "integral")
```

gets the same adjusted p-values and p-value matrix from
`adjust_iwt(t0, t_coeff, recycle = FALSE)`.

### fd objects

Wherever fdatest takes curves, fdatestcpp takes an `fd` object or a matrix
(one curve per row, one grid point per column), and handles `fd` objects as
fdatest does:

- an `fd` object is evaluated on `seq(range[1], range[2], by = dx)` with
  `t(fda::eval.fd(...))`. The default `dx = NULL` uses 1% of the basis range,
  i.e. 101 points. `dx` is ignored for matrices;
- in two-sample tests, both samples must be `fd` objects on the same range,
  or both matrices;
- `mu` can be a number, a vector with one value per grid point or, with `fd`
  data, an `fd` object holding one function on the same range;
- the output's `data` field holds the evaluated curves, and an `fd` `mu` is
  returned as the `1 x p` matrix of its values, as in fdatest. The extra
  field `abscissa` holds the grid.

## Exactness

fdatest computes its statistics with R primitives (`colMeans()`, `rowSums()`,
`sum()`, `cov()`, `mean()`) that accumulate in 80-bit long double and round to
double at the end. fdatestcpp accumulates in `long double` in the same order
and rounds at the same points, including `sum()`'s overflow to `Inf`. It draws
signs as `unif_rand() >= 0.5`, which is what `rbinom(1, 1, 0.5)` returns from
the same uniform, and permutations with `R_unif_index()` as `sample()` does,
in fdatest's order, so the random number stream ends in the same state.
p-values are `count / B` where fdatest divides counts, and long double means
where it calls `mean()`. The IWT computes every interval statistic by running
sums from each start point, which reproduces R's rounding (differences of
prefix sums would not); the TWT sums a level set in grid order whenever a
faster sorted-order sum is too close to the observed value to decide the
comparison, or close to overflow.

As a result, with the same inputs and seed, unadjusted p-values, adjusted
p-values and IWT p-value matrices are `identical()` to fdatest's, for matrices
and `fd` objects, with one deliberate exception: a two-sample `mu` that varies
along the grid, which fdatest centres incorrectly (see
[Differences from fdatest](#differences-from-fdatest)).

This needs a compiler whose `long double` has the precision of R's (as with
the usual x86_64 builds of R on Linux, Windows and Intel macOS);
`fdatestcpp:::cpp_acc_digits()` and `.Machine$longdouble.digits` show both.
It has been verified on Windows x86_64 (see [Tests](#tests)). Where
`long double` is plain double, as on macOS with Apple silicon, R's own sums
change in the same way, so exact agreement is likely there too but has not
been verified.

### Tests

The scripts in `tests/` run with `R CMD check`, or one by one with
`Rscript tests/test-engine.R` from the package folder after installing the
package. They check the IWT adjustment against a brute-force maximum over the
p-value matrix, the TWT against a direct grid-order computation, the
resampling against plain R computations, `fd` input against its
evaluation passed as a matrix, and the plots (`test-plot.R`): the shaded
cells against the adjusted p-values, and every kind of result, panel and
option built without warnings (685 checks).

The comparisons with fdatest itself are opt-in, so that a future fdatest
release cannot break this package's checks. Set the environment variable
`FDATESTCPP_COMPARE_FDATEST=true` to run them; they need the development
version of fdatest and the long double condition above. They require
`identical()` results for:

- the engine (`test-engine.R`), on statistics from real tests, uniform draws,
  values with ties, `NaN` and `Inf` entries, p = 1 and B = 1, for both
  aggregations and both values of `recycle`;
- the resampling (`test-resample.R`): `t0`, `t_perm` and the final random
  number state against fdatest's one- and two-sample permutation code, for
  all alternatives, paired and unpaired, with and without `standardize`;
- the wrappers (`test-wrappers.R`) and `fd` input (`test-fd.R`: Fourier and
  B-spline bases, ranges not starting at 0, default and explicit `dx`,
  scalar, vector and `fd` `mu`);
- edge cases found while reviewing the package (`test-review.R`): names,
  overflow near the largest double, nearly constant curves with
  `standardize`, user-supplied permutations, a single block.

With R 4.6.1 on Windows x86_64, ggplot2 4.0.3, fdatest at commit `67bb41c`
and fda 6.3.0, the six scripts run 1,980 checks with the comparisons,
with no failures. They also pass when the package is compiled with the
undefined-behaviour sanitizer and `_GLIBCXX_ASSERTIONS`, and with `double`
instead of `long double` accumulators (the comparisons are then skipped).

## Speed

- **IWT.** For each resample, every interval statistic is accumulated from
  its start point one grid point at a time, so the p(p+1)/2 intervals (p^2
  with `recycle = TRUE`) cost O(p^2) instead of fdatest's O(p^3) R loops. The
  adjusted p-values come from an O(p^2) recursion.
- **TWT.** Sorting the grid points by unadjusted p-value turns the level sets
  into prefixes and suffixes, so all set statistics of a resample come out of
  one pass: O(B p) instead of O(B p^2).
- **Resampling.** One pass over the curves per resample; with `blocks = q`
  the curves are first summed within blocks, and each resample costs O(q p).
- **fd input** adds one evaluation of the curves on the grid.

The code is single-threaded: to use several cores, parallelize across
simulation replications or datasets.

End-to-end times for n = 240 curves (or pairs), p = 200 grid points,
B = 5000 resamples and `recycle = FALSE`, on one core of an Intel Core
i7-1360P with R 4.6.1:

| Test | fdatest, integral | fdatestcpp, integral | fdatest, max | fdatestcpp, max |
|---|---:|---:|---:|---:|
| `iwt1` | 35 s | 0.22 s (162x) | 153 s | 0.14 s (1,074x) |
| `iwt2(paired = TRUE)` | 34 s | 0.29 s (117x) | 147 s | 0.25 s (588x) |
| `twt2(paired = TRUE)` | 4.4 s | 0.22 s (20x) | 6.2 s | 0.21 s (29x) |

Over a grid of 36 settings (p from 40 to 200, B of 1000 and 5000) the
speed-up ranges from 17x to 1,074x, with median 61x, and every setting gave
`identical()` p-values. With `fd` input (240 curves on fdatest's default
101-point grid, B of 1000 and 5000) the speed-up ranges from 18x to 469x.
The full tables and the scripts that produce them are in
[benchmarks.md](benchmarks.md) and `bench/`.

## Differences from fdatest

- **Two-sample `mu` that varies along the grid** (a vector, or a
  non-constant `fd` function): fdatest subtracts `matrix(mu, n1, p)`, which
  fills by column, so its centering is wrong unless `mu` is constant.
  fdatestcpp subtracts `mu[j]` at grid point `j`. Results agree when `mu`
  takes exactly the same value at every grid point.
- **Stricter input checks.** Data must be numeric, finite and non-empty
  (fdatest also takes logical matrices and drops `NA`s in two-sample means);
  paired samples must have the same number of curves; flags such as `paired`
  and `recycle` must be `TRUE` or `FALSE`; `mu` must be finite with length 1
  or p; an `fd` `mu` must hold one function. Inputs fdatest fails on with an
  unclear error (matrix data with an `fd` `mu`, multivariate `fd` objects)
  stop with a message naming the argument.
- **Inputs accepted that fdatest rejects:** data frames, `mu` as a `1 x p`
  matrix (such as a previous result's `mu`), abbreviated `alternative` and
  `aggregation_strategy`.
- **Unpaired tests on a one-point grid** work (fdatest's fail because
  `coeff[sample(n), ]` drops the matrix dimension).
- **No parallel resampling.** fdatest can draw its resamples in parallel with
  mirai daemons; its results then depend on the daemons and cannot be
  matched seed for seed.
- **Global p-value of the IWT.** It is `pvalue_matrix[1, 1]` (count / B, as in
  fdatest's matrix), while `global1()` and `global2()` use `mean()`, as
  fdatest's `global2()`. The two can differ in the last bit when B >= 2048,
  and with `NaN` statistics the IWT ignores them while the global test
  returns `NA`.
- **Plots.** fdatestcpp has its own `plot()` and `autoplot()` methods (see
  [Plots](#plots)); fdatest's methods remain available as
  `fdatest:::plot.fos()` and `fdatest:::plot.fts()`.
- **Not ported:** ANOVA and functional linear model tests, PCT, FDR and the
  `ITP*` functions, and the `IWTimage()` heatmap (fdatest's works on the
  outputs).

## References

Pini, A. and Vantini, S. (2017). Interval-wise testing for functional data.
*Journal of Nonparametric Statistics*, 29(2), 407-424.
[doi:10.1080/10485252.2017.1306627](https://doi.org/10.1080/10485252.2017.1306627)

Abramowicz, K., Pini, A., Schelin, L., Sjöstedt de Luna, S., Stamm, A. and
Vantini, S. (2023). Domain selection and familywise error rate for functional
data: a unified framework. *Biometrics*, 79(2), 1119-1132.
[doi:10.1111/biom.13669](https://doi.org/10.1111/biom.13669)

Pini, A., Stamm, A. and Vantini, S. fdatest: Interval Wise Testing for
Functional Data. R package, development version (commit `67bb41c`).
<https://github.com/permaverse/fdatest>

## Authors and license

fdatestcpp is written and maintained by Matteo Fontana. The testing
procedures, and the fdatest code that fdatestcpp ports, are by Alessia Pini,
Aymeric Stamm and Simone Vantini, who are listed as contributors and
copyright holders. Please cite the papers above when you use the tests.

License: GPL (>= 3), as the development version of fdatest.
