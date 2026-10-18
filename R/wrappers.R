# Drop-in versions of fdatest's iwt1(), iwt2(), twt2() and global2() (same
# argument names, defaults and output fields; matrices or fd objects), plus
# one-sample twt1() and global1(). Extra arguments: alternative/standardize for
# the one-sample tests, blocks, signs, perms, plus_one and return_perm. The
# outputs have fdatest's classes "fos"/"fts", preceded by fdatestcpp's own
# subclass for the plot methods in plot.R.

new_fos <- function(prep, adj, correction, st, return_perm) {
  out <- list(data = prep$data, mu = prep$mu,
              unadjusted_pvalues = adj$unadjusted_pvalues,
              adjusted_pvalues = adj$adjusted_pvalues,
              correction_method = correction)
  if (correction == "IWT") out$pvalue_matrix <- adj$pvalue_matrix
  out$global_pvalue <- adj$global_pvalue
  out$t0 <- st$t0
  if (check_flag(return_perm)) out$t_perm <- st$t_perm
  if (!is.null(prep$abscissa)) out$abscissa <- prep$abscissa
  class(out) <- c("fdatestcpp_fos", "fos")
  out
}

new_fts <- function(prep, adj, correction, st, return_perm) {
  out <- list(data = prep$data,
              group_labels = c(rep(1, nrow(prep$x1)), rep(2, nrow(prep$x2))),
              mu = prep$mu,
              unadjusted_pvalues = adj$unadjusted_pvalues,
              adjusted_pvalues = adj$adjusted_pvalues,
              correction_method = correction)
  if (correction == "IWT") out$pvalue_matrix <- adj$pvalue_matrix
  out$global_pvalue <- adj$global_pvalue
  out$t0 <- st$t0
  if (check_flag(return_perm)) out$t_perm <- st$t_perm
  if (!is.null(prep$abscissa)) out$abscissa <- prep$abscissa
  class(out) <- c("fdatestcpp_fts", "fts")
  out
}

iwt1 <- function(data, mu = 0, n_perm = 1000L, dx = NULL, recycle = TRUE,
                 verbose = FALSE, aggregation_strategy = c("integral", "max"),
                 alternative = c("two.sided", "less", "greater"), standardize = FALSE,
                 blocks = NULL, signs = NULL, plus_one = FALSE, return_perm = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  prep <- prep_one(data, mu, dx)
  st <- stats_one(prep, n_perm, match.arg(alternative), standardize, blocks, signs)
  adj <- adjust_iwt(st$t0, st$t_perm, recycle, aggregation_strategy, plus_one)
  new_fos(prep, adj, "IWT", st, return_perm)
}

twt1 <- function(data, mu = 0, n_perm = 1000L, dx = NULL, verbose = FALSE,
                 aggregation_strategy = c("integral", "max"),
                 alternative = c("two.sided", "less", "greater"), standardize = FALSE,
                 blocks = NULL, signs = NULL, plus_one = FALSE, return_perm = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  prep <- prep_one(data, mu, dx)
  st <- stats_one(prep, n_perm, match.arg(alternative), standardize, blocks, signs)
  adj <- adjust_twt(st$t0, st$t_perm, aggregation_strategy, plus_one)
  new_fos(prep, adj, "TWT", st, return_perm)
}

global1 <- function(data, mu = 0, n_perm = 1000L, dx = NULL, verbose = FALSE,
                    aggregation_strategy = c("integral", "max"),
                    alternative = c("two.sided", "less", "greater"), standardize = FALSE,
                    blocks = NULL, signs = NULL, plus_one = FALSE, return_perm = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  prep <- prep_one(data, mu, dx)
  st <- stats_one(prep, n_perm, match.arg(alternative), standardize, blocks, signs)
  adj <- adjust_global(st$t0, st$t_perm, aggregation_strategy, plus_one)
  new_fos(prep, adj, "Global", st, return_perm)
}

iwt2 <- function(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
                 alternative = c("two.sided", "less", "greater"), standardize = FALSE,
                 verbose = FALSE, aggregation_strategy = c("integral", "max"),
                 recycle = TRUE, blocks = NULL, signs = NULL, perms = NULL,
                 plus_one = FALSE, return_perm = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  prep <- prep_two(data1, data2, mu, dx)
  st <- stats_two(prep, n_perm, paired, match.arg(alternative), standardize, blocks,
                  signs, perms)
  adj <- adjust_iwt(st$t0, st$t_perm, recycle, aggregation_strategy, plus_one)
  new_fts(prep, adj, "IWT", st, return_perm)
}

twt2 <- function(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
                 alternative = c("two.sided", "less", "greater"), standardize = FALSE,
                 verbose = FALSE, aggregation_strategy = c("integral", "max"),
                 blocks = NULL, signs = NULL, perms = NULL,
                 plus_one = FALSE, return_perm = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  prep <- prep_two(data1, data2, mu, dx)
  st <- stats_two(prep, n_perm, paired, match.arg(alternative), standardize, blocks,
                  signs, perms)
  adj <- adjust_twt(st$t0, st$t_perm, aggregation_strategy, plus_one)
  new_fts(prep, adj, "TWT", st, return_perm)
}

global2 <- function(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
                    alternative = c("two.sided", "less", "greater"), standardize = FALSE,
                    verbose = FALSE, aggregation_strategy = c("integral", "max"),
                    blocks = NULL, signs = NULL, perms = NULL,
                    plus_one = FALSE, return_perm = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  prep <- prep_two(data1, data2, mu, dx)
  st <- stats_two(prep, n_perm, paired, match.arg(alternative), standardize, blocks,
                  signs, perms)
  adj <- adjust_global(st$t0, st$t_perm, aggregation_strategy, plus_one)
  new_fts(prep, adj, "Global", st, return_perm)
}
