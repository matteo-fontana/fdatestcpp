# Timing: fdatest vs fdatestcpp, end to end, on the same data and seed.
#
# n = 240 curves (or pairs), p in {40, 100, 200} grid points,
# B in {1000, 5000} resamples; iwt1, iwt2(paired = TRUE), twt2(paired = TRUE);
# aggregation "integral" and "max"; recycle = FALSE (non-periodic domain).
# fdatest runs once per cell (median of three when a run takes under 5 s);
# fdatestcpp is timed as the median of repeated calls. Pin the process to the
# fast cores on hybrid CPUs (Windows: start /affinity FF Rscript ...). Each cell also checks that the two packages return identical
# unadjusted and adjusted p-values (and IWT p-value matrices).
#
# Run from the fdatestcpp folder:  Rscript bench/bench_fdatest_vs_port.R
# Writes bench/results.csv (appended cell by cell; finished cells are skipped
# on a rerun) and ../benchmarks.md via bench/write_benchmarks_md.R.

suppressPackageStartupMessages({
  library(data.table)
  library(fdatest)
  library(fdatestcpp)
})
options(warn = -1)  # fdatest warns on max() over empty TWT sets
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script)) setwd(dirname(dirname(normalizePath(script))))

out_csv <- file.path("bench", "results.csv")
n <- 240L

# Simulated curves: FAR(1) noise with an exponential
# kernel, a bump alternative on tau > 0.6, and a persistent common level.
make_data <- function(p, psi = 0.5, seed = 1L) {
  set.seed(seed)
  tau <- (seq_len(p) - 0.5) / p
  R <- chol(exp(-abs(outer(tau, tau, "-")) / 0.3))
  far1 <- function(psi, burn = 100L) {
    Z <- matrix(rnorm((n + burn) * p), n + burn, p) %*% R
    E <- Z
    for (t in 2:(n + burn)) E[t, ] <- psi * E[t - 1, ] + sqrt(1 - psi^2) * Z[t, ]
    E[(burn + 1):(burn + n), , drop = FALSE]
  }
  mu <- ifelse(tau > 0.6, 3.5 * (0.5 + tau) * sqrt((1 + psi) / (1 - psi)) / sqrt(n), 0)
  d <- sweep(sweep(far1(psi), 2, 0.5 + tau, "*"), 2, mu, "+")
  common <- exp(0.5 * far1(0.9))
  list(d = d, l1 = common + d / 2, l2 = common - d / 2)
}

call_test <- function(pkg, fun, agg, dat, B) {
  ns <- asNamespace(pkg)
  switch(fun,
    iwt1 = ns$iwt1(dat$d, mu = 0, n_perm = B, recycle = FALSE, aggregation_strategy = agg),
    iwt2_paired = ns$iwt2(dat$l1, dat$l2, n_perm = B, paired = TRUE, recycle = FALSE,
                          aggregation_strategy = agg),
    twt2_paired = ns$twt2(dat$l1, dat$l2, n_perm = B, paired = TRUE,
                          aggregation_strategy = agg))
}

# fdatest: one run, or the median of three when a run takes less than 5 s
time_fdatest <- function(f) {
  gc(FALSE)
  t <- system.time(res <- f())[["elapsed"]]
  if (t < 5) {
    t <- median(c(t, vapply(1:2, function(i) {
      gc(FALSE)
      system.time(f())[["elapsed"]]
    }, numeric(1))))
  }
  list(time = t, result = res)
}

# median seconds per call; fast calls are looped until a batch lasts >= 0.25 s
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
  }, numeric(1)))
}

grid <- CJ(p = c(40L, 100L, 200L), B = c(1000L, 5000L),
           fun = c("iwt1", "iwt2_paired", "twt2_paired"), agg = c("integral", "max"),
           sorted = FALSE)
done <- if (file.exists(out_csv)) fread(out_csv) else NULL

# warm-up: load and byte-compile everything once
warm <- make_data(10L)
for (pkg in c("fdatest", "fdatestcpp")) for (f in unique(grid$fun)) for (a in unique(grid$agg)) {
  invisible(call_test(pkg, f, a, warm, 20L))
}

for (i in seq_len(nrow(grid))) {
  g <- grid[i]
  if (!is.null(done) && nrow(done[p == g$p & B == g$B & fun == g$fun & agg == g$agg])) next
  dat <- make_data(g$p)
  seed <- 1000L + i
  fdt <- time_fdatest(function() {
    set.seed(seed)
    call_test("fdatest", g$fun, g$agg, dat, g$B)
  })
  t_fd <- fdt$time
  r_fd <- fdt$result
  set.seed(seed)
  r_cpp <- call_test("fdatestcpp", g$fun, g$agg, dat, g$B)
  same <- identical(r_fd$unadjusted_pvalues, r_cpp$unadjusted_pvalues) &&
    identical(r_fd$adjusted_pvalues, r_cpp$adjusted_pvalues) &&
    identical(r_fd$pvalue_matrix, r_cpp$pvalue_matrix)
  t_cpp <- time_port(function() {
    set.seed(seed)
    call_test("fdatestcpp", g$fun, g$agg, dat, g$B)
  })
  row <- data.table(p = g$p, B = g$B, fun = g$fun, agg = g$agg, n = n,
                    fdatest_s = t_fd, fdatestcpp_s = t_cpp, speedup = t_fd / t_cpp,
                    identical = same, n_reject_adj = sum(r_cpp$adjusted_pvalues <= 0.05))
  fwrite(row, out_csv, append = file.exists(out_csv))
  cat(sprintf("p=%3d B=%4d %-12s %-8s fdatest %8.2fs  fdatestcpp %7.4fs  x%7.0f  identical=%s\n",
              g$p, g$B, g$fun, g$agg, t_fd, t_cpp, t_fd / t_cpp, same))
}

source(file.path("bench", "write_benchmarks_md.R"))
