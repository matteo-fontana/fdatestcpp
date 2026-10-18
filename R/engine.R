# P-value engine on precomputed pointwise statistics.

# fdatest's unadjusted p-values, colSums(t_perm >= t0) / B, carry the column
# names of t_perm; its TWT adjusted p-values start from them and keep them.
name_pvalues <- function(out, t_perm, twt = FALSE) {
  nm <- colnames(t_perm)
  if (!is.null(nm)) {
    names(out$unadjusted_pvalues) <- nm
    if (twt) names(out$adjusted_pvalues) <- nm
  }
  out
}

pvalues_pointwise <- function(t0, t_perm, plus_one = FALSE) {
  s <- check_stats(t0, t_perm)
  out <- cpp_pointwise(s$t0, s$t_perm, check_flag(plus_one))
  names(out) <- colnames(t_perm)
  out
}

adjust_iwt <- function(t0, t_perm, recycle, aggregation_strategy = c("integral", "max"),
                       plus_one = FALSE, pvalue_matrix = TRUE) {
  if (missing(recycle)) {
    stop(paste0("`recycle` must be given: TRUE for a periodic domain, FALSE otherwise ",
                "(fdatest's iwt1()/iwt2() default to TRUE)."), call. = FALSE)
  }
  aggregation_strategy <- match.arg(aggregation_strategy)
  s <- check_stats(t0, t_perm)
  out <- cpp_iwt(s$t0, s$t_perm, check_flag(recycle), aggregation_strategy == "max",
                 check_flag(plus_one), check_flag(pvalue_matrix))
  name_pvalues(out, t_perm)
}

adjust_twt <- function(t0, t_perm, aggregation_strategy = c("integral", "max"),
                       plus_one = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  s <- check_stats(t0, t_perm)
  out <- cpp_twt(s$t0, s$t_perm, aggregation_strategy == "max", check_flag(plus_one), FALSE)
  name_pvalues(out, t_perm, twt = TRUE)
}

adjust_global <- function(t0, t_perm, aggregation_strategy = c("integral", "max"),
                          plus_one = FALSE) {
  aggregation_strategy <- match.arg(aggregation_strategy)
  s <- check_stats(t0, t_perm)
  out <- cpp_global(s$t0, s$t_perm, aggregation_strategy == "max", check_flag(plus_one))
  name_pvalues(out, t_perm)
}
