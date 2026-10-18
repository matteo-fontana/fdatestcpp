# Timing of fdatestcpp's stages at n = 240, p = 200, B = 5000.
# Run from the fdatestcpp folder:  Rscript bench/bench_stages.R
# Writes bench/stages.csv; bench/write_benchmarks_md.R adds it to benchmarks.md.

suppressPackageStartupMessages({
  library(data.table)
  library(fdatestcpp)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script)) setwd(dirname(dirname(normalizePath(script))))
n <- 240L
p <- 200L
B <- 5000L
set.seed(1)
tau <- (seq_len(p) - 0.5) / p
R <- chol(exp(-abs(outer(tau, tau, "-")) / 0.3))
l1 <- matrix(rnorm(n * p), n, p) %*% R
l2 <- matrix(rnorm(n * p), n, p) %*% R
d <- l1 - l2

# median seconds per call; fast calls are looped until a batch lasts >= 0.25 s
med_time <- function(f, batches = 5L) {
  k <- 1L
  repeat {
    t <- system.time(for (i in seq_len(k)) f())[["elapsed"]]
    if (t >= 0.25 || k >= 1024L) break
    k <- k * 4L
  }
  median(vapply(seq_len(batches), function(i) {
    gc(FALSE)
    system.time(for (i in seq_len(k)) f())[["elapsed"]] / k
  }, numeric(1)))
}

s1 <- perm_stats1(d, n_perm = B)
stages <- list(
  "Resampling: one-sample sign flips" = function() perm_stats1(d, n_perm = B),
  "Resampling: one-sample, 12 blocks" = function() perm_stats1(d, n_perm = B, blocks = 12),
  "Resampling: paired swaps" = function() perm_stats2(l1, l2, n_perm = B, paired = TRUE),
  "Resampling: paired swaps, standardized" =
    function() perm_stats2(l1, l2, n_perm = B, paired = TRUE, standardize = TRUE),
  "Resampling: unpaired permutations" = function() perm_stats2(l1, l2, n_perm = B),
  "IWT, integral, recycle = FALSE" = function() adjust_iwt(s1$t0, s1$t_perm, FALSE, "integral"),
  "IWT, max, recycle = FALSE" = function() adjust_iwt(s1$t0, s1$t_perm, FALSE, "max"),
  "IWT, integral, recycle = TRUE" = function() adjust_iwt(s1$t0, s1$t_perm, TRUE, "integral"),
  "TWT, integral" = function() adjust_twt(s1$t0, s1$t_perm, "integral"),
  "TWT, max" = function() adjust_twt(s1$t0, s1$t_perm, "max"),
  "Global, integral" = function() adjust_global(s1$t0, s1$t_perm, "integral")
)
res <- data.table(stage = names(stages), seconds = vapply(stages, med_time, numeric(1)))
fwrite(res, file.path("bench", "stages.csv"))
print(res)
