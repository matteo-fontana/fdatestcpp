# Regression tests for problems found while reviewing the fd support
# (October 2026): names of p-values, dx with matrix data, perms
# validation, overflow near DBL_MAX, and the one-sample studentized statistic.

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
  all(vapply(c("iwt1", "iwt2", "twt2", "global2", "os_permtest", "p_adjust_iwt", "p_adjust_twt", "p_adjust_global"), exists, NA,
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
pv <- c("unadjusted_pvalues", "adjusted_pvalues")

set.seed(7)
n <- 15L
p <- 12L
X1 <- matrix(rnorm(n * p), n, p) + 0.3
X2 <- matrix(rnorm(n * p), n, p)
X3 <- matrix(rnorm((n + 4L) * p), n + 4L, p)
B <- 150L

# 1. Column names: fdatest's colMeans()/colSums() name the unadjusted p-values
#    (and the TWT adjusted ones); rbind() takes the first matrix's names
if (has_fdatest) {
  nmA <- paste0("q", seq_len(p))
  nmB <- paste0("m", seq_len(p))
  A1 <- X1; colnames(A1) <- nmA
  A2 <- X2; colnames(A2) <- nmB
  A3 <- X3; colnames(A3) <- nmB
  for (agg in c("integral", "max")) {
    set.seed(1); a <- fdatest::iwt1(A1, n_perm = B, recycle = FALSE, aggregation_strategy = agg)
    set.seed(1); b <- fdatestcpp::iwt1(A1, n_perm = B, recycle = FALSE, aggregation_strategy = agg)
    same(a, b, paste("colnames iwt1", agg), c(pv, "pvalue_matrix", "data"))
    for (cs in list(list(A1, A2, "both"), list(X1, A2, "data2 only"), list(A1, X2, "data1 only"))) {
      for (paired in c(TRUE, FALSE)) {
        y <- if (paired) cs[[2]] else if (is.null(colnames(cs[[2]]))) X3 else A3
        lab <- sprintf("colnames %s paired=%s %s", cs[[3]], paired, agg)
        set.seed(2); a <- fdatest::iwt2(cs[[1]], y, n_perm = B, paired = paired,
                                        aggregation_strategy = agg, recycle = FALSE)
        set.seed(2); b <- fdatestcpp::iwt2(cs[[1]], y, n_perm = B, paired = paired,
                                           aggregation_strategy = agg, recycle = FALSE)
        same(a, b, paste("iwt2", lab), c(pv, "pvalue_matrix", "data"))
        set.seed(3); a <- fdatest::twt2(cs[[1]], y, n_perm = B, paired = paired,
                                        aggregation_strategy = agg, standardize = TRUE)
        set.seed(3); b <- fdatestcpp::twt2(cs[[1]], y, n_perm = B, paired = paired,
                                           aggregation_strategy = agg, standardize = TRUE)
        same(a, b, paste("twt2", lab), pv)
        set.seed(4); a <- fdatest::global2(cs[[1]], y, n_perm = B, paired = paired,
                                           aggregation_strategy = agg, alternative = "less")
        set.seed(4); b <- fdatestcpp::global2(cs[[1]], y, n_perm = B, paired = paired,
                                              aggregation_strategy = agg, alternative = "less")
        same(a, b, paste("global2", lab), pv)
      }
    }
  }
  # internals: t0 and t_perm carry the names as in os_permtest()/ts_permtest()
  set.seed(5); a <- asNamespace("fdatest")$os_permtest(A1, B)
  set.seed(5); b <- perm_stats1(A1, n_perm = B)
  expect(identical(a$t0, b$t0) && identical(a$t_coeff, b$t_perm), "colnames os_permtest")

  # 2. dx is ignored for matrix data, whatever its value
  for (dx in list(-1, 0, NA, "a", c(0.1, 0.2))) {
    set.seed(6); a <- fdatest::iwt1(X1, n_perm = B, dx = dx, recycle = FALSE)
    set.seed(6); b <- fdatestcpp::iwt1(X1, n_perm = B, dx = dx, recycle = FALSE)
    same(a, b, sprintf("matrix data, dx = %s", format(dx)[1]), c(pv, "pvalue_matrix"))
  }
}
expect(!err(twt2(X1, X2, dx = 0, n_perm = 20)), "dx ignored for matrices, two-sample")

# 3. perms: whole numbers in 1..N only, names ignored, no coercion warnings
P <- t(replicate(10, sample(2L * n + 4L)))
expect(err(perm_stats2(X1, X3, perms = P + 0.7)), "non-integer perms rejected")
expect(!err(perm_stats2(X1, X3, perms = P)), "integer-valued perms accepted")
Pn <- P
colnames(Pn) <- paste0("obs", seq_len(ncol(P)))
expect(identical(perm_stats2(X1, X3, perms = Pn), perm_stats2(X1, X3, perms = P)),
       "perms with column names")
for (bad in c(Inf, -Inf, 1e10, 0, 2L * n + 5L)) {
  Pb <- P
  Pb[1, 1] <- bad
  w <- 0L
  e <- withCallingHandlers(err(perm_stats2(X1, X3, perms = Pb)),
                           warning = function(cond) { w <<- w + 1L; invokeRestart("muffleWarning") })
  expect(e && w == 0L, sprintf("perms entry %s rejected without warnings", format(bad)))
}

# 3b. blocks = 1 (one block) equals the id vector rep(1, n)
S1 <- matrix(sample(c(-1, 1), 20, TRUE), 20, 1)
expect(identical(perm_stats1(X1, signs = S1, blocks = 1),
                 perm_stats1(X1, signs = S1, blocks = rep(1, n))), "blocks = 1")

# 4. Overflow near DBL_MAX: fdatest sums T0 with sum() (beyond DBL_MAX -> Inf)
#    and T* with rowSums() (rounded); TWT decides overflowing sets exactly
if (has_fdatest) {
  fd <- asNamespace("fdatest")
  M <- .Machine$double.xmax
  cases <- list(
    list(t0 = rep(1e308, 4), tp = rbind(c(0, 0, 1.7e308, 1.7e308), c(0, 0, 0, 0),
                                        c(0, 0, 0, 1.7e308), c(0, 0, 0, 0)), lab = "sums overflow"),
    list(t0 = c(M, 2^969), tp = rbind(c(M, 2^968), c(0, 0)), lab = "DBL_MAX window"),
    list(t0 = -c(M, 2^969), tp = rbind(-c(M, 2^968), c(0, 0), -c(M, 0)), lab = "negative window"))
  set.seed(8)
  big <- matrix(runif(60 * 8, 0.9, 1) * 2^1021, 60, 8)
  cases[[4]] <- list(t0 = big[1, ], tp = big[-1, ], lab = "random near overflow")
  for (cs in cases) {
    t0 <- cs$t0; tp <- cs$tp; Bc <- nrow(tp); pc <- ncol(tp)
    pval <- colSums(tp >= matrix(t0, Bc, pc, byrow = TRUE)) / Bc
    for (rec in c(FALSE, TRUE)) {
      a <- fd$p_adjust_iwt(pc, pval, t0, tp, Bc, rec, "integral")
      b <- adjust_iwt(t0, tp, rec)
      expect(identical(a$adjusted_pvalues, b$adjusted_pvalues) &&
               identical(a$pvalue_matrix, b$pvalue_matrix), paste(cs$lab, "IWT", rec))
    }
    a <- fd$p_adjust_twt(pval, pc, t0, tp, "integral")
    expect(identical(a$adjusted_pvalues, adjust_twt(t0, tp)$adjusted_pvalues), paste(cs$lab, "TWT"))
    a <- fd$p_adjust_global("integral", t0, tp, pc)
    expect(identical(a$adjusted_pvalues, adjust_global(t0, tp)$adjusted_pvalues),
           paste(cs$lab, "Global"))
  }
  expect(identical(adjust_twt(cases[[1]]$t0, cases[[1]]$tp)$global_pvalue,
                   adjust_global(cases[[1]]$t0, cases[[1]]$tp)$global_pvalue),
         "TWT and Global agree on the full set under overflow")
}

# 5. One-sample standardize = mean^2 / var() of the flipped curves, exactly as
#    R computes colMeans(y)^2 / apply(y, 2, var), also where curves are equal,
#    nearly equal or nearly symmetric (a one-pass variance fails there)
# element-wise agreement (identical on x86_64, where both use 80-bit long
# double; within rounding elsewhere), with the same Inf and NaN positions
close <- function(a, b, tol = 1e-6) {
  fin <- is.finite(a) & is.finite(b)
  identical(is.na(a), is.na(b)) && identical(is.infinite(a), is.infinite(b)) &&
    identical(a[is.infinite(a)], b[is.infinite(b)]) &&
    all(abs(a[fin] - b[fin]) <= tol * pmax(abs(b[fin]), .Machine$double.xmin))
}
ref_std <- function(x, S, id, alt = "two.sided") {
  t(apply(S, 1, function(s) {
    y <- x * s[id]
    m <- colMeans(y)
    num <- switch(alt, two.sided = m^2, greater = (m * (m > 0))^2, less = (m * (m < 0))^2)
    num / apply(y, 2, var)
  }))
}
set.seed(9)
nn <- 12L
Z <- matrix(rnorm(nn * 6), nn, 6)
Z[, 1] <- 0.3                                        # all curves equal
Z[, 2] <- 2.2 + 1e-12 * rnorm(nn)                    # nearly equal
Z[, 3] <- 1e8 + rnorm(nn)                            # large level, unit noise
Z[, 4] <- 4.2e7 + rnorm(nn)
Z[, 5] <- c(1.34573317902666, -1.3457331790244966, rep(0, nn - 2))  # near symmetry
for (blk in list(NULL, 4L)) {
  q <- if (is.null(blk)) nn else blk
  id <- if (is.null(blk)) seq_len(nn) else cut(seq_len(nn), blk, labels = FALSE)
  S <- matrix(sample(c(-1, 1), 200 * q, TRUE), 200, q)
  S[1, ] <- 1
  S[2, ] <- -1
  for (alt in c("two.sided", "greater", "less")) {
    b <- perm_stats1(Z, signs = S, blocks = blk, standardize = TRUE, alternative = alt)
    r <- ref_std(Z, S, id, alt)
    # columns 1-5 are ill-conditioned: compare them with R only when both
    # accumulate with the same precision; column 6 is ordinary data
    cols <- if (exact_platform) seq_len(ncol(Z)) else 6L
    expect(close(b$t_perm[, cols], r[, cols]),
           sprintf("standardized flips = R var(), blocks %s, %s", format(blk), alt))
    expect(close(b$t0[cols], r[1, cols]),
           sprintf("standardized t0, blocks %s, %s", format(blk), alt))
    expect(identical(is.na(b$t_perm), is.na(r)) &&
             identical(is.infinite(b$t_perm), is.infinite(r)),
           sprintf("standardized Inf/NaN positions, blocks %s, %s", format(blk), alt))
    expect(all(b$t_perm[!is.na(b$t_perm)] >= 0), sprintf("no negative statistics, %s", alt))
  }
}
b <- perm_stats1(Z, signs = matrix(1, 1, nn), standardize = TRUE)
expect(is.infinite(b$t0[1]) && b$t0[1] > 0, "equal curves: t0 = +Inf")
expect(all(b$t0[-6] > 0), "no negative studentized statistics")

# 6. p = 1 unpaired (fdatest crashes; the port works)
expect(!err(iwt2(X1[, 1, drop = FALSE], X3[, 1, drop = FALSE], n_perm = 20)), "p = 1 unpaired")

cat(sprintf("test-review: %d checks, %d failed%s\n", n_checks, length(fails),
            if (has_fdatest) "" else " (fdatest comparisons skipped)"))
if (length(fails)) stop("Failed:\n", paste(fails, collapse = "\n"))
