# Timing with fd input: fdatest vs fdatestcpp, end to end, same data and seed.
#
# n = 240 fd curves (pairs for iwt2/twt2) on a 20-function B-spline basis,
# evaluated on fdatest's default grid (dx = 1% of the range, 101 points);
# B in {1000, 5000}; iwt1, iwt2(paired = TRUE), twt2(paired = TRUE); both
# aggregations; recycle = FALSE. Also the cost of the fd evaluation itself
# (fd input vs the same curves passed as a matrix).
#
# Run from the fdatestcpp folder, pinned to fast cores on hybrid CPUs:
#   start "" /affinity FF /wait /b Rscript bench\bench_fd.R
# Writes bench/fd_results.csv; bench/write_benchmarks_md.R adds it to benchmarks.md.

suppressPackageStartupMessages({
  library(data.table)
  library(fda)
  library(fdatest)
  library(fdatestcpp)
})
options(warn = -1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script)) setwd(dirname(dirname(normalizePath(script))))

n <- 240L
# simulated curves (as in bench_fdatest_vs_port.R) on 60 points, smoothed
set.seed(1)
tt <- (seq_len(60) - 0.5) / 60
R <- chol(exp(-abs(outer(tt, tt, "-")) / 0.3))
far1 <- function(psi, burn = 100L) {
  Z <- matrix(rnorm((n + burn) * 60), n + burn, 60) %*% R
  E <- Z
  for (t in 2:(n + burn)) E[t, ] <- psi * E[t - 1, ] + sqrt(1 - psi^2) * Z[t, ]
  E[(burn + 1):(burn + n), , drop = FALSE]
}
mu <- ifelse(tt > 0.6, 3.5 * (0.5 + tt) * sqrt(3) / sqrt(n), 0)
d <- sweep(sweep(far1(0.5), 2, 0.5 + tt, "*"), 2, mu, "+")
common <- exp(0.5 * far1(0.9))
basis <- create.bspline.basis(c(0, 1), 20)
smooth <- function(y) smooth.basis(tt, t(y), basis)$fd
fdat <- list(d = smooth(d), l1 = smooth(common + d / 2), l2 = smooth(common - d / 2))
grid <- seq(0, 1, by = 0.01)
mdat <- lapply(fdat, function(f) t(eval.fd(grid, f)))

call_test <- function(pkg, fun, agg, dat, B) {
  ns <- asNamespace(pkg)
  switch(fun,
    iwt1 = ns$iwt1(dat$d, mu = 0, n_perm = B, recycle = FALSE, aggregation_strategy = agg),
    iwt2_paired = ns$iwt2(dat$l1, dat$l2, n_perm = B, paired = TRUE, recycle = FALSE,
                          aggregation_strategy = agg),
    twt2_paired = ns$twt2(dat$l1, dat$l2, n_perm = B, paired = TRUE,
                          aggregation_strategy = agg))
}
time_once_or_3 <- function(f) {
  gc(FALSE)
  t <- system.time(res <- f())[["elapsed"]]
  if (t < 5) t <- median(c(t, vapply(1:2, function(i) system.time(f())[["elapsed"]], 1)))
  list(time = t, result = res)
}
time_port <- function(f, batches = 5L) {
  k <- 1L
  repeat {
    t <- system.time(for (i in seq_len(k)) f())[["elapsed"]]
    if (t >= 0.25 || k >= 256L) break
    k <- k * 4L
  }
  median(vapply(seq_len(batches), function(i) {
    gc(FALSE)
    system.time(for (i in seq_len(k)) f())[["elapsed"]] / k
  }, 1))
}

for (pkg in c("fdatest", "fdatestcpp")) for (f in c("iwt1", "iwt2_paired", "twt2_paired")) {
  invisible(call_test(pkg, f, "integral", fdat, 20L))
}
grid_cells <- CJ(B = c(1000L, 5000L), fun = c("iwt1", "iwt2_paired", "twt2_paired"),
                 agg = c("integral", "max"), sorted = FALSE)
rows <- list()
for (i in seq_len(nrow(grid_cells))) {
  g <- grid_cells[i]
  seed <- 2000L + i
  fdt <- time_once_or_3(function() { set.seed(seed); call_test("fdatest", g$fun, g$agg, fdat, g$B) })
  set.seed(seed)
  r_cpp <- call_test("fdatestcpp", g$fun, g$agg, fdat, g$B)
  same <- identical(fdt$result$unadjusted_pvalues, r_cpp$unadjusted_pvalues) &&
    identical(fdt$result$adjusted_pvalues, r_cpp$adjusted_pvalues) &&
    identical(fdt$result$pvalue_matrix, r_cpp$pvalue_matrix) &&
    identical(fdt$result$data, r_cpp$data)
  t_cpp <- time_port(function() { set.seed(seed); call_test("fdatestcpp", g$fun, g$agg, fdat, g$B) })
  t_mat <- time_port(function() { set.seed(seed); call_test("fdatestcpp", g$fun, g$agg, mdat, g$B) })
  rows[[i]] <- data.table(p = length(grid), B = g$B, fun = g$fun, agg = g$agg, n = n,
                          fdatest_s = fdt$time, fdatestcpp_s = t_cpp, matrix_input_s = t_mat,
                          speedup = fdt$time / t_cpp, identical = same)
  cat(sprintf("fd B=%4d %-12s %-8s fdatest %7.2fs  fdatestcpp %7.4fs (matrix %7.4fs)  x%5.0f  identical=%s\n",
              g$B, g$fun, g$agg, fdt$time, t_cpp, t_mat, fdt$time / t_cpp, same))
}
fwrite(rbindlist(rows), file.path("bench", "fd_results.csv"))
source(file.path("bench", "write_benchmarks_md.R"))
