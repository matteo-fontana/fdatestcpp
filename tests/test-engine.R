# Engine on precomputed statistics.
# Always: internal consistency (sorted-order TWT = index-order TWT, IWT
# correction = max over intervals containing each point, pointwise p-values).
# With the development version of fdatest (which has p_adjust_iwt; the CRAN
# version does not): every output identical() to fdatest:::p_adjust_iwt /
# p_adjust_twt / p_adjust_global on the same t0, t_coeff and pval.

library(fdatestcpp)
# Exact comparisons with fdatest are opt-in (FDATESTCPP_COMPARE_FDATEST=true),
# so that changes in fdatest cannot break this package's checks. They need
# fdatest's development version (the CRAN release 2.1.1 lacks these objects)
# and an accumulator with R's long double precision (not so, e.g., on R
# built without long double), see fdatestcpp:::cpp_acc_digits().
r_digits <- if (capabilities("long.double")) .Machine$longdouble.digits else 53L
exact_platform <- identical(as.integer(r_digits), fdatestcpp:::cpp_acc_digits())
has_fdatest <- identical(Sys.getenv("FDATESTCPP_COMPARE_FDATEST"), "true") &&
  exact_platform && requireNamespace("fdatest", quietly = TRUE) &&
  all(vapply(c("p_adjust_iwt", "p_adjust_twt", "p_adjust_global", "pval_correct_cpp"), exists, NA,
             envir = asNamespace("fdatest"), inherits = FALSE))
fd <- if (has_fdatest) asNamespace("fdatest") else NULL
if (has_fdatest) options(warn = -1)  # fdatest warns on max() over empty TWT sets

n_checks <- 0L
fails <- character(0)
expect <- function(ok, label) {
  n_checks <<- n_checks + 1L
  if (!isTRUE(ok)) fails <<- c(fails, label)
}

# Pointwise statistics of several kinds
gen <- function(kind, p, B) {
  if (kind == "flip") {                       # sign-flip test, bump alternative
    n <- 30L
    tau <- seq_len(p) / p
    x <- matrix(rnorm(n * p), n, p) %*% chol(exp(-abs(outer(tau, tau, "-")) / 0.3))
    x <- sweep(x, 2, ifelse(tau > 0.6, 0.4, 0), "+")
    t0 <- colMeans(x)^2
    tp <- t(replicate(B, colMeans(x * sample(c(-1, 1), n, TRUE))^2))
    tp[sample(B, 3), ] <- matrix(t0, 3, p, byrow = TRUE)   # identity flips: exact ties
  } else if (kind == "unif") {
    t0 <- runif(p)
    tp <- matrix(runif(B * p), B, p)
  } else if (kind == "discrete") {            # ties in value, summation order matters
    v <- c(0, 0.1, 0.2, 0.3)
    t0 <- sample(v, p, TRUE)
    tp <- matrix(sample(v, B * p, TRUE), B, p)
  } else if (kind == "signed") {              # negative statistics (user-supplied)
    t0 <- rnorm(p)
    tp <- matrix(rnorm(B * p), B, p)
    tp[1, ] <- t0
  } else if (kind == "nonfinite") {
    t0 <- runif(p)
    tp <- matrix(runif(B * p), B, p)
    tp[2, 3] <- NaN
    tp[4, p] <- Inf
    tp[5, 2] <- Inf
    tp[5, 4] <- -Inf
    t0[p - 1] <- NaN
  }
  list(t0 = t0, tp = tp)
}

# IWT correction by brute force: max (ignoring NA) of the p-values of all
# intervals containing each point; row i of the matrix holds intervals of
# length p - i + 1 starting at column j
naive_correct <- function(pm, recycle) {
  p <- nrow(pm)
  adj <- rep(NA_real_, p)
  for (i in seq_len(p)) {
    len <- p - i + 1L
    for (j in if (recycle) seq_len(p) else seq_len(i)) {
      v <- pm[i, j]
      if (is.na(v)) next
      pts <- ((j - 1L + seq_len(len) - 1L) %% p) + 1L
      adj[pts] <- ifelse(is.na(adj[pts]) | adj[pts] < v, v, adj[pts])
    }
  }
  adj
}

compare_all <- function(s, label) {
  t0 <- s$t0
  tp <- s$tp
  B <- nrow(tp)
  p <- ncol(tp)
  pv <- colSums(tp >= matrix(t0, B, p, byrow = TRUE)) / B   # as fdatest's permtest
  expect(identical(pvalues_pointwise(t0, tp), pv), paste(label, "pointwise"))
  for (agg in c("integral", "max")) {
    for (rec in c(FALSE, TRUE)) {
      b <- adjust_iwt(t0, tp, rec, agg)
      lab <- sprintf("%s IWT %s recycle=%s", label, agg, rec)
      expect(identical(b$adjusted_pvalues, naive_correct(b$pvalue_matrix, rec)),
             paste(lab, "correction = brute force"))
      expect(identical(b$pvalue_matrix[p, ], b$unadjusted_pvalues), paste(lab, "last row"))
    }
    b <- adjust_twt(t0, tp, agg)
    bb <- fdatestcpp:::cpp_twt(t0, tp, agg == "max", FALSE, TRUE)   # index-order sums
    expect(identical(b$adjusted_pvalues, bb$adjusted_pvalues),
           sprintf("%s TWT %s sorted = index order", label, agg))
  }
  if (has_fdatest) compare_fdatest(s, label, pv)
}

compare_fdatest <- function(s, label, pv) {
  t0 <- s$t0
  tp <- s$tp
  B <- nrow(tp)
  p <- ncol(tp)
  for (agg in c("integral", "max")) {
    for (rec in c(FALSE, TRUE)) {
      a <- fd$p_adjust_iwt(p, pv, t0, tp, B, rec, agg)
      b <- adjust_iwt(t0, tp, rec, agg)
      lab <- sprintf("%s IWT %s recycle=%s", label, agg, rec)
      expect(identical(b$unadjusted_pvalues, pv), paste(lab, "unadjusted"))
      expect(identical(b$adjusted_pvalues, a$adjusted_pvalues), paste(lab, "adjusted"))
      expect(identical(b$pvalue_matrix, a$pvalue_matrix), paste(lab, "pvalue_matrix"))
      expect(identical(b$global_pvalue, a$pvalue_matrix[1, 1]), paste(lab, "global"))
      expect(identical(fdatestcpp:::cpp_iwt_correct(a$pvalue_matrix, rec),
                       fd$pval_correct_cpp(a$pvalue_matrix)[1, ]), paste(lab, "correction"))
    }
    a <- fd$p_adjust_twt(pv, p, t0, tp, agg)
    b <- adjust_twt(t0, tp, agg)
    lab <- sprintf("%s TWT %s", label, agg)
    expect(identical(b$unadjusted_pvalues, pv), paste(lab, "unadjusted"))
    expect(identical(b$adjusted_pvalues, a$adjusted_pvalues), paste(lab, "adjusted"))
    bb <- fdatestcpp:::cpp_twt(t0, tp, agg == "max", FALSE, TRUE)   # index-order sums
    expect(identical(bb$adjusted_pvalues, a$adjusted_pvalues), paste(lab, "adjusted, brute"))
    a <- fd$p_adjust_global(agg, t0, tp, p)
    b <- adjust_global(t0, tp, agg)
    lab <- sprintf("%s Global %s", label, agg)
    expect(identical(b$adjusted_pvalues, a$adjusted_pvalues), paste(lab, "adjusted"))
  }
}

set.seed(20261016)
for (kind in c("flip", "unif", "discrete", "signed")) {
  for (p in c(7L, 40L, 101L)) {
    B <- if (p == 101L) 100L else 200L
    compare_all(gen(kind, p, B), sprintf("[%s p=%d B=%d]", kind, p, B))
  }
}
compare_all(gen("nonfinite", 12L, 80L), "[nonfinite p=12]")
compare_all(gen("unif", 1L, 50L), "[p=1]")
compare_all(gen("unif", 9L, 1L), "[B=1]")
s <- gen("unif", 10L, 60L)
compare_all(list(t0 = s$t0, tp = matrix(s$t0, 60, 10, byrow = TRUE)), "[all tied]")
compare_all(list(t0 = rep(0, 10), tp = matrix(0, 60, 10)), "[all zero]")

# plus_one: (1 + count) / (B + 1)
s <- gen("flip", 30L, 200L)
a <- adjust_iwt(s$t0, s$tp, FALSE)
b <- adjust_iwt(s$t0, s$tp, FALSE, plus_one = TRUE)
expect(isTRUE(all.equal(b$pvalue_matrix, (a$pvalue_matrix * 200 + 1) / 201)), "plus_one IWT matrix")
expect(isTRUE(all.equal(b$adjusted_pvalues, (a$adjusted_pvalues * 200 + 1) / 201)), "plus_one IWT adjusted")
a <- adjust_twt(s$t0, s$tp)
b <- adjust_twt(s$t0, s$tp, plus_one = TRUE)
expect(isTRUE(all.equal(b$adjusted_pvalues, (a$adjusted_pvalues * 200 + 1) / 201)), "plus_one TWT")
a <- adjust_global(s$t0, s$tp)
b <- adjust_global(s$t0, s$tp, plus_one = TRUE)
expect(isTRUE(all.equal(b$global_pvalue, (a$global_pvalue * 200 + 1) / 201)), "plus_one Global")

# argument checks
expect(inherits(try(adjust_iwt(s$t0, s$tp), silent = TRUE), "try-error"), "recycle is required")
expect(inherits(try(adjust_iwt(s$t0[-1], s$tp, FALSE), silent = TRUE), "try-error"), "dimension check")

cat(sprintf("test-engine: %d checks, %d failed%s\n", n_checks, length(fails),
            if (has_fdatest) "" else " (fdatest comparisons skipped)"))
if (length(fails)) stop("Failed:\n", paste(fails, collapse = "\n"))
