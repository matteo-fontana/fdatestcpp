# Front end: observed and resampled pointwise statistics.

perm_stats1 <- function(data, mu = 0, n_perm = 1000L, dx = NULL,
                        alternative = c("two.sided", "less", "greater"),
                        standardize = FALSE, blocks = NULL, signs = NULL) {
  prep <- prep_one(data, mu, dx)
  out <- stats_one(prep, n_perm, match.arg(alternative), standardize, blocks, signs)
  out$mu <- prep$mu
  out$abscissa <- prep$abscissa
  out
}

perm_stats2 <- function(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
                        alternative = c("two.sided", "less", "greater"),
                        standardize = FALSE, blocks = NULL, signs = NULL, perms = NULL) {
  prep <- prep_two(data1, data2, mu, dx)
  out <- stats_two(prep, n_perm, paired, match.arg(alternative), standardize, blocks,
                   signs, perms)
  out$mu <- prep$mu
  out$abscissa <- prep$abscissa
  out
}

# Statistics from prepared data (see prepare.R)
stats_one <- function(prep, n_perm, alternative, standardize, blocks, signs) {
  x <- prep$x
  n <- nrow(x)
  p <- ncol(x)
  coeff <- x - matrix(prep$mu, nrow = n, ncol = p, byrow = TRUE)
  units <- make_units(blocks, n)
  if (is.null(signs)) {
    B <- check_nperm(n_perm)
    ind <- matrix(0L, 0L, units$q)
  } else {
    signs <- check_signs(signs, units$q)
    B <- nrow(signs)
    ind <- matrix(as.integer(signs > 0), B, units$q)
  }
  out <- cpp_stats_onesample(coeff, units$id, units$q, B, ind, alt_code(alternative),
                             check_flag(standardize))
  out$n_units <- units$q
  name_stats(out, colnames(coeff))
}

# Grid-point names as fdatest's colMeans() gives them: the column names of the
# (first named) data matrix, if any; the p-values inherit them in the engine.
name_stats <- function(out, nm) {
  if (!is.null(nm)) {
    names(out$t0) <- nm
    colnames(out$t_perm) <- nm
  }
  out
}

stats_two <- function(prep, n_perm, paired, alternative, standardize, blocks, signs,
                      perms) {
  data1 <- prep$x1
  data2 <- prep$x2
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p <- ncol(data1)
  # mu is subtracted from column j of data1. fdatest uses matrix(mu, n1, p),
  # which fills by column: the same for a constant mu, not for a varying one.
  x1 <- data1 - matrix(prep$mu, nrow = n1, ncol = p, byrow = TRUE)
  alt <- alt_code(alternative)
  standardize <- check_flag(standardize)
  if (check_flag(paired)) {
    if (n1 != n2) {
      stop("A paired test needs `data1` and `data2` with the same number of curves.",
           call. = FALSE)
    }
    if (!is.null(perms)) {
      stop("`perms` is for unpaired tests; paired swaps take `signs` (-1 = swap).",
           call. = FALSE)
    }
    units <- make_units(blocks, n1)
    if (is.null(signs)) {
      B <- check_nperm(n_perm)
      ind <- matrix(0L, 0L, units$q)
    } else {
      signs <- check_signs(signs, units$q)
      B <- nrow(signs)
      ind <- matrix(as.integer(signs < 0), B, units$q)
    }
    out <- if (standardize) {
      cpp_stats_twosample_obs(rbind(x1, data2), n1, TRUE, units$id, units$q, B, ind,
                              alt, TRUE)
    } else {
      cpp_stats_paired(x1, data2, units$id, units$q, B, ind, alt)
    }
    out$n_units <- units$q
  } else {
    if (!is.null(blocks)) {
      stop("`blocks` applies to one-sample sign flips and paired swaps only.", call. = FALSE)
    }
    if (!is.null(signs)) {
      stop("`signs` applies to sign flips and paired swaps; unpaired tests take `perms`.",
           call. = FALSE)
    }
    N <- n1 + n2
    if (is.null(perms)) {
      B <- check_nperm(n_perm)
      ind <- matrix(0L, 0L, N)
    } else {
      ind <- check_perms(perms, N)
      B <- nrow(ind)
    }
    out <- cpp_stats_twosample_obs(rbind(x1, data2), n1, FALSE, integer(0), 0L, B, ind,
                                   alt, standardize)
    out$n_units <- N
  }
  # rbind(data1, data2) takes the column names of the first matrix that has them
  name_stats(out, if (is.null(colnames(x1))) colnames(data2) else colnames(x1))
}
