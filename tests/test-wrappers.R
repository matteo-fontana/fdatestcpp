# Wrappers end to end. Always: wrappers = front end + engine on the same draws.
# With the development version of fdatest (which has iwt1; the CRAN version
# does not): with the same seed, fdatestcpp::iwt1/iwt2/twt2/global2 return
# identical() p-values and p-value matrices to fdatest's.

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

set.seed(2)
n <- 30L
p <- 24L
B <- 250L
tau <- (seq_len(p) - 0.5) / p
R <- chol(exp(-abs(outer(tau, tau, "-")) / 0.3))
d <- matrix(rnorm(n * p), n, p) %*% R + matrix(ifelse(tau > 0.6, 0.5, 0), n, p, byrow = TRUE)
common <- matrix(rnorm(n * p), n, p) %*% R
l1 <- common + d / 2
l2 <- common - d / 2
other <- matrix(rnorm((n + 5L) * p), n + 5L, p) %*% R

pv_fields <- c("unadjusted_pvalues", "adjusted_pvalues")

# wrappers = front end + engine, with the same draws
for (agg in c("integral", "max")) {
  set.seed(9); st <- perm_stats1(d, mu = 0.05, n_perm = B)
  set.seed(9); b <- iwt1(d, mu = 0.05, n_perm = B, recycle = FALSE, aggregation_strategy = agg)
  e <- adjust_iwt(st$t0, st$t_perm, FALSE, agg)
  same(b, e, paste("iwt1 = perm_stats1 + adjust_iwt", agg), c(pv_fields, "pvalue_matrix"))
  for (paired in c(TRUE, FALSE)) {
    y <- if (paired) l2 else other
    set.seed(9); st <- perm_stats2(l1, y, n_perm = B, paired = paired, standardize = TRUE)
    set.seed(9); b <- twt2(l1, y, n_perm = B, paired = paired, standardize = TRUE,
                           aggregation_strategy = agg)
    same(b, adjust_twt(st$t0, st$t_perm, agg), paste("twt2 = perm_stats2 + adjust_twt", agg),
         pv_fields)
  }
}

if (has_fdatest) for (agg in c("integral", "max")) {
  for (rec in c(FALSE, TRUE)) {
    set.seed(10); a <- fdatest::iwt1(d, mu = 0.05, n_perm = B, recycle = rec,
                                     aggregation_strategy = agg)
    set.seed(10); b <- fdatestcpp::iwt1(d, mu = 0.05, n_perm = B, recycle = rec,
                                        aggregation_strategy = agg)
    same(a, b, sprintf("iwt1 %s recycle=%s", agg, rec),
         c(pv_fields, "pvalue_matrix", "data", "mu"))
    expect(identical(class(b), c("fdatestcpp_fos", class(a))), "iwt1 class")
    for (paired in c(TRUE, FALSE)) {
      y <- if (paired) l2 else other
      for (alt in c("two.sided", "greater", "less")) {
        set.seed(11); a <- fdatest::iwt2(l1, y, n_perm = B, paired = paired, alternative = alt,
                                         aggregation_strategy = agg, recycle = rec)
        set.seed(11); b <- fdatestcpp::iwt2(l1, y, n_perm = B, paired = paired,
                                            alternative = alt, aggregation_strategy = agg,
                                            recycle = rec)
        same(a, b, sprintf("iwt2 paired=%s %s %s recycle=%s", paired, alt, agg, rec),
             c(pv_fields, "pvalue_matrix", "data", "group_labels", "mu"))
      }
    }
  }
  for (paired in c(TRUE, FALSE)) {
    y <- if (paired) l2 else other
    for (std in c(FALSE, TRUE)) {
      set.seed(12); a <- fdatest::twt2(l1, y, n_perm = B, paired = paired, standardize = std,
                                       aggregation_strategy = agg)
      set.seed(12); b <- fdatestcpp::twt2(l1, y, n_perm = B, paired = paired,
                                          standardize = std, aggregation_strategy = agg)
      same(a, b, sprintf("twt2 paired=%s standardize=%s %s", paired, std, agg), pv_fields)
      set.seed(13); a <- fdatest::global2(l1, y, n_perm = B, paired = paired,
                                          standardize = std, aggregation_strategy = agg)
      set.seed(13); b <- fdatestcpp::global2(l1, y, n_perm = B, paired = paired,
                                             standardize = std, aggregation_strategy = agg)
      same(a, b, sprintf("global2 paired=%s standardize=%s %s", paired, std, agg),
           c(pv_fields, "global_pvalue"))
      set.seed(14); a <- fdatest::iwt2(l1, y, n_perm = B, paired = paired, standardize = std,
                                       aggregation_strategy = agg, recycle = FALSE)
      set.seed(14); b <- fdatestcpp::iwt2(l1, y, n_perm = B, paired = paired,
                                          standardize = std, aggregation_strategy = agg,
                                          recycle = FALSE)
      same(a, b, sprintf("iwt2 paired=%s standardize=%s %s", paired, std, agg),
           c(pv_fields, "pvalue_matrix"))
    }
  }
}

# One-sample TWT / Global (not in fdatest) agree with the engine on the same draws
set.seed(15); st <- perm_stats1(d, n_perm = B)
set.seed(15); b <- twt1(d, n_perm = B, return_perm = TRUE)
expect(identical(b$adjusted_pvalues, adjust_twt(st$t0, st$t_perm)$adjusted_pvalues), "twt1")
expect(identical(b$t_perm, st$t_perm), "twt1 return_perm")
set.seed(15); b <- global1(d, n_perm = B)
expect(identical(b$global_pvalue, adjust_global(st$t0, st$t_perm)$global_pvalue), "global1")

cat(sprintf("test-wrappers: %d checks, %d failed%s\n", n_checks, length(fails),
            if (has_fdatest) "" else " (fdatest comparisons skipped)"))
if (length(fails)) stop("Failed:\n", paste(fails, collapse = "\n"))
