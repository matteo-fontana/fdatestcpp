# Benchmarks: fdatest vs fdatestcpp

Run on 2026-10-16: R version 4.6.1 (2026-06-24 ucrt), Windows 11, Intel Core i7-1360P laptop (Intel64 Family 6 Model 186 Stepping 2, GenuineIntel), single thread pinned to the performance cores (`start /affinity FF`), no mirai daemons.
fdatest 2.1.1 from https://alessiapini.r-universe.dev, commit 67bb41c; fdatestcpp 0.1.0.

End to end, same data and seed: n = 240 curves (pairs for iwt2/twt2), p grid points, B resamples, `recycle = FALSE`; simulated curves (FAR(1) noise, a bump in the mean on part of the domain). fdatest is timed once per cell, or as the median of three runs when a run takes under 5 s; fdatestcpp is the median over repeated calls. "Identical" means `identical()` unadjusted and adjusted p-values and, for the IWT, the whole p-value matrix.

| Test | Aggregation | p | B | fdatest (s) | fdatestcpp (s) | Speed-up | Identical p-values |
|:---|:---|---:|---:|---:|---:|---:|:---|
| iwt1 | integral | 40 | 1000 | 0.110 | 0.0064 | 17x | yes |
| iwt1 | integral | 40 | 5000 | 0.580 | 0.027 | 22x | yes |
| iwt1 | integral | 100 | 1000 | 0.920 | 0.014 | 64x | yes |
| iwt1 | integral | 100 | 5000 | 5.0 | 0.075 | 67x | yes |
| iwt1 | integral | 200 | 1000 | 7.2 | 0.044 | 165x | yes |
| iwt1 | integral | 200 | 5000 | 35 | 0.218 | 162x | yes |
| iwt1 | max | 40 | 1000 | 0.730 | 0.0055 | 133x | yes |
| iwt1 | max | 40 | 5000 | 4.1 | 0.024 | 169x | yes |
| iwt1 | max | 100 | 1000 | 6.3 | 0.015 | 429x | yes |
| iwt1 | max | 100 | 5000 | 33 | 0.068 | 487x | yes |
| iwt1 | max | 200 | 1000 | 31 | 0.034 | 897x | yes |
| iwt1 | max | 200 | 5000 | 153 | 0.142 | 1,074x | yes |
| iwt2 (paired) | integral | 40 | 1000 | 0.330 | 0.010 | 32x | yes |
| iwt2 (paired) | integral | 40 | 5000 | 1.7 | 0.058 | 29x | yes |
| iwt2 (paired) | integral | 100 | 1000 | 1.5 | 0.026 | 58x | yes |
| iwt2 (paired) | integral | 100 | 5000 | 6.6 | 0.130 | 51x | yes |
| iwt2 (paired) | integral | 200 | 1000 | 7.6 | 0.060 | 126x | yes |
| iwt2 (paired) | integral | 200 | 5000 | 34 | 0.290 | 117x | yes |
| iwt2 (paired) | max | 40 | 1000 | 1.1 | 0.0098 | 112x | yes |
| iwt2 (paired) | max | 40 | 5000 | 5.3 | 0.047 | 113x | yes |
| iwt2 (paired) | max | 100 | 1000 | 6.9 | 0.024 | 284x | yes |
| iwt2 (paired) | max | 100 | 5000 | 33 | 0.132 | 252x | yes |
| iwt2 (paired) | max | 200 | 1000 | 32 | 0.055 | 576x | yes |
| iwt2 (paired) | max | 200 | 5000 | 147 | 0.250 | 588x | yes |
| twt2 (paired) | integral | 40 | 1000 | 0.320 | 0.011 | 29x | yes |
| twt2 (paired) | integral | 40 | 5000 | 1.4 | 0.046 | 31x | yes |
| twt2 (paired) | integral | 100 | 1000 | 0.570 | 0.023 | 25x | yes |
| twt2 (paired) | integral | 100 | 5000 | 2.5 | 0.112 | 22x | yes |
| twt2 (paired) | integral | 200 | 1000 | 0.900 | 0.046 | 20x | yes |
| twt2 (paired) | integral | 200 | 5000 | 4.4 | 0.220 | 20x | yes |
| twt2 (paired) | max | 40 | 1000 | 0.380 | 0.012 | 33x | yes |
| twt2 (paired) | max | 40 | 5000 | 2.0 | 0.049 | 41x | yes |
| twt2 (paired) | max | 100 | 1000 | 0.790 | 0.026 | 30x | yes |
| twt2 (paired) | max | 100 | 5000 | 3.6 | 0.110 | 32x | yes |
| twt2 (paired) | max | 200 | 1000 | 1.5 | 0.052 | 28x | yes |
| twt2 (paired) | max | 200 | 5000 | 6.2 | 0.212 | 29x | yes |

Speed-up range: 17x to 1,074x; median 61x.
Total time for the grid: fdatest 579 s, fdatestcpp 2.7 s.

The IWT gains grow with p because fdatest's adjustment loops over all p(p+1)/2 intervals in R (`rowSums()` on a copied submatrix for the integral, `apply(..., 1, max)` for the max), O(B p^3) work with a large constant. The TWT gains are smaller because fdatest's TWT adjustment is already vectorized over resamples; there the port's time is mostly the paired resampling (see below), which fdatest does with an R loop over pairs.

## fd input

n = 240 fd curves (B-spline basis with 20 functions) evaluated on fdatest's default grid (101 points), `recycle = FALSE`, same seed. The "same curves as matrix" column times fdatestcpp on the evaluated curves, so the difference to the fd column is the cost of the evaluation. "Identical" also covers the evaluated `data` field.

| Test | Aggregation | B | fdatest (s) | fdatestcpp (s) | same curves as matrix (s) | Speed-up | Identical p-values |
|:---|:---|---:|---:|---:|---:|---:|:---|
| iwt1 | integral | 1000 | 0.840 | 0.017 | 0.015 | 50x | yes |
| iwt1 | integral | 5000 | 4.4 | 0.070 | 0.067 | 62x | yes |
| iwt1 | max | 1000 | 5.3 | 0.012 | 0.012 | 425x | yes |
| iwt1 | max | 5000 | 28 | 0.060 | 0.058 | 469x | yes |
| iwt2 (paired) | integral | 1000 | 1.3 | 0.026 | 0.026 | 49x | yes |
| iwt2 (paired) | integral | 5000 | 5.9 | 0.130 | 0.125 | 46x | yes |
| iwt2 (paired) | max | 1000 | 6.0 | 0.024 | 0.024 | 244x | yes |
| iwt2 (paired) | max | 5000 | 30 | 0.115 | 0.115 | 257x | yes |
| twt2 (paired) | integral | 1000 | 0.430 | 0.024 | 0.022 | 18x | yes |
| twt2 (paired) | integral | 5000 | 2.3 | 0.110 | 0.110 | 21x | yes |
| twt2 (paired) | max | 1000 | 0.630 | 0.022 | 0.021 | 28x | yes |
| twt2 (paired) | max | 5000 | 3.2 | 0.105 | 0.105 | 31x | yes |

## Where fdatestcpp spends its time

Single calls, n = 240, p = 200, B = 5000 (median of repeated calls). Resampling is the front end (`perm_stats1()` / `perm_stats2()`); the adjustments run on its output.

| Stage | Seconds |
|:---|---:|
| Resampling: one-sample sign flips | 0.090 |
| Resampling: one-sample, 12 blocks | 0.0089 |
| Resampling: paired swaps | 0.222 |
| Resampling: paired swaps, standardized | 1.3 |
| Resampling: unpaired permutations | 0.420 |
| IWT, integral, recycle = FALSE | 0.087 |
| IWT, max, recycle = FALSE | 0.047 |
| IWT, integral, recycle = TRUE | 0.182 |
| TWT, integral | 0.0083 |
| TWT, max | 0.0052 |
| Global, integral | 0.0038 |
