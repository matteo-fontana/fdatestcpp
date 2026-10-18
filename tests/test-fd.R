# fd objects: fdatestcpp evaluates them as fdatest does. Always: fd input
# equals its evaluation passed as a matrix, and argument errors. With the
# development version of fdatest: with the same seed the wrappers return
# identical() results to fdatest's for fd input.

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
  all(vapply(c("iwt1", "iwt2", "twt2", "global2"), exists, NA,
             envir = asNamespace("fdatest"), inherits = FALSE))
if (has_fdatest) options(warn = -1)  # fdatest warns on max() over empty TWT sets

n_checks <- 0L
fails <- character(0)
expect <- function(ok, label) {
  n_checks <<- n_checks + 1L
  if (!isTRUE(ok)) fails <<- c(fails, label)
}
same <- function(a, b, label, fields) {
  for (f in fields) expect(identical(a[[f]], b[[f]]), paste(label, f))
}
err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

# fd samples: smoothed noisy curves on a Fourier basis (periodic, [0, 12]) and
# on a B-spline basis (non-periodic, [-2, 3.5])
set.seed(3)
make_fd <- function(n, basis, shift = function(t) 0, sd = 1) {
  rv <- basis$rangeval
  tt <- seq(rv[1], rv[2], length.out = 60)
  y <- sapply(seq_len(n), function(i) {
    shift(tt) + cumsum(rnorm(60, sd = sd)) / 4 + rnorm(1)
  })
  fda::smooth.basis(tt, y, basis)$fd
}
fb <- fda::create.fourier.basis(c(0, 12), 13)
bb <- fda::create.bspline.basis(c(-2, 3.5), 12)
F1 <- make_fd(20, fb, function(t) ifelse(t > 8, 1, 0))
F2 <- make_fd(20, fb)
F3 <- make_fd(26, fb)
S1 <- make_fd(18, bb, function(t) 0.8 * (t > 1))
S2 <- make_fd(18, bb)
mu_f <- fda::fd(matrix(c(0.6, 0.4, -0.2), 3, 1), fda::create.fourier.basis(c(0, 12), 3))
mu_const <- fda::fd(matrix(c(0.5, 0, 0), 3, 1), fda::create.fourier.basis(c(0, 12), 3))
B <- 120L
pv <- c("unadjusted_pvalues", "adjusted_pvalues")

if (has_fdatest) {
  cases1 <- list(
    list(x = F1, mu = 0, dx = NULL, lab = "fourier default dx"),
    list(x = F1, mu = 0.3, dx = 0.5, lab = "fourier dx 0.5"),
    list(x = F1, mu = sin(seq(0, 12, by = 0.5)), dx = 0.5, lab = "fourier vector mu"),
    list(x = F1, mu = mu_f, dx = NULL, lab = "fourier fd mu"),
    list(x = S1, mu = 0, dx = 0.07, lab = "bspline dx 0.07"),
    list(x = S1, mu = 0.1, dx = NULL, lab = "bspline default dx"))
  for (cs in cases1) {
    for (agg in c("integral", "max")) for (rec in c(FALSE, TRUE)) {
      set.seed(21); a <- fdatest::iwt1(cs$x, mu = cs$mu, n_perm = B, dx = cs$dx,
                                       recycle = rec, aggregation_strategy = agg)
      set.seed(21); b <- fdatestcpp::iwt1(cs$x, mu = cs$mu, n_perm = B, dx = cs$dx,
                                          recycle = rec, aggregation_strategy = agg)
      same(a, b, sprintf("iwt1 %s %s recycle=%s", cs$lab, agg, rec),
           c(pv, "pvalue_matrix", "data", "mu", "correction_method"))
      expect(identical(class(b), c("fdatestcpp_fos", class(a))), "iwt1 class")
    }
  }
  cases2 <- list(
    list(x = F1, y = F2, mu = 0, dx = NULL, lab = "fourier"),
    list(x = F1, y = F2, mu = rep(0.2, 41), dx = 0.3, lab = "fourier constant vector mu"),
    list(x = F1, y = F2, mu = mu_const, dx = 0.3, lab = "fourier const fd mu dx 0.3"),
    list(x = S1, y = S2, mu = 0.2, dx = 0.11, lab = "bspline dx 0.11"))
  for (cs in cases2) {
    for (paired in c(TRUE, FALSE)) {
      y <- if (paired || !identical(cs$y, F2)) cs$y else F3
      for (agg in c("integral", "max")) {
        lab <- sprintf("%s paired=%s %s", cs$lab, paired, agg)
        set.seed(22); a <- fdatest::iwt2(cs$x, y, mu = cs$mu, dx = cs$dx, n_perm = B,
                                         paired = paired, aggregation_strategy = agg,
                                         recycle = FALSE)
        set.seed(22); b <- fdatestcpp::iwt2(cs$x, y, mu = cs$mu, dx = cs$dx, n_perm = B,
                                            paired = paired, aggregation_strategy = agg,
                                            recycle = FALSE)
        same(a, b, paste("iwt2", lab), c(pv, "pvalue_matrix", "data", "mu", "group_labels"))
        set.seed(23); a <- fdatest::twt2(cs$x, y, mu = cs$mu, dx = cs$dx, n_perm = B,
                                         paired = paired, aggregation_strategy = agg,
                                         standardize = TRUE, alternative = "greater")
        set.seed(23); b <- fdatestcpp::twt2(cs$x, y, mu = cs$mu, dx = cs$dx, n_perm = B,
                                            paired = paired, aggregation_strategy = agg,
                                            standardize = TRUE, alternative = "greater")
        same(a, b, paste("twt2 standardized greater", lab), c(pv, "data", "mu"))
        set.seed(24); a <- fdatest::global2(cs$x, y, mu = cs$mu, dx = cs$dx, n_perm = B,
                                            paired = paired, aggregation_strategy = agg)
        set.seed(24); b <- fdatestcpp::global2(cs$x, y, mu = cs$mu, dx = cs$dx, n_perm = B,
                                               paired = paired, aggregation_strategy = agg)
        same(a, b, paste("global2", lab), c(pv, "global_pvalue", "data", "mu"))
      }
    }
  }
}

# fd input = its evaluation on fdatest's grid, passed as a matrix
grid <- seq(0, 12, by = 0.12)
X1 <- t(fda::eval.fd(grid, F1))
X2 <- t(fda::eval.fd(grid, F2))
set.seed(31); a <- iwt2(F1, F2, paired = TRUE, n_perm = B, recycle = FALSE, blocks = 5)
set.seed(31); b <- iwt2(X1, X2, paired = TRUE, n_perm = B, recycle = FALSE, blocks = 5)
same(a, b, "fd = evaluated matrix (paired, blocks)", c(pv, "pvalue_matrix", "data"))
expect(identical(a$abscissa, grid), "abscissa of the default grid")
expect(is.null(b$abscissa), "no abscissa for matrix input")
s <- perm_stats1(F1, dx = 0.5, n_perm = 10)
expect(identical(s$abscissa, seq(0, 12, by = 0.5)) && ncol(s$t_perm) == 25L,
       "perm_stats1 grid with dx")
# non-constant fd mu in two-sample tests: mu(t_j) is subtracted at grid point j
mu_vals <- t(fda::eval.fd(grid, mu_f))
set.seed(32); a <- twt2(F1, F2, mu = mu_f, paired = TRUE, n_perm = B)
set.seed(32); b <- twt2(X1 - matrix(mu_vals, nrow(X1), ncol(X1), byrow = TRUE), X2,
                        paired = TRUE, n_perm = B)
same(a, b, "non-constant fd mu, two-sample", pv)
expect(identical(a$mu, mu_vals), "fd mu returned as its values")

# errors
expect(err(iwt2(F1, X2, n_perm = 10)), "fd + matrix")
expect(err(iwt1(X1, mu = mu_f, n_perm = 10)), "matrix data + fd mu")
expect(err(iwt2(F1, S2, n_perm = 10)), "different ranges")
expect(err(iwt1(S1, mu = mu_f, n_perm = 10)), "mu range differs")
expect(err(iwt1(F1, dx = -1, n_perm = 10)), "negative dx")
msg <- tryCatch(iwt1(F1, dx = 1e-12, n_perm = 10), error = conditionMessage)
expect(is.character(msg) && grepl("dx", msg), "dx too small: message names dx")
two_mu <- fda::fd(matrix(c(1, 0, 0, 2, 0, 0), 3, 2), fda::create.fourier.basis(c(0, 12), 3))
expect(err(iwt1(F1, mu = two_mu, n_perm = 10)), "mu with two curves")
mv <- fda::fd(array(rnorm(13 * 4 * 2), c(13, 4, 2)), fb)
expect(err(iwt1(mv, n_perm = 10)), "multivariate fd")

cat(sprintf("test-fd: %d checks, %d failed%s\n", n_checks, length(fails),
            if (has_fdatest) "" else " (fdatest comparisons skipped)"))
if (length(fails)) stop("Failed:\n", paste(fails, collapse = "\n"))
