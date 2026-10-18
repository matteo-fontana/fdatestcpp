# Argument checks shared by the engine, the front end and the wrappers.

check_flag <- function(x, name = deparse(substitute(x))) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

check_nperm <- function(n_perm) {
  if (!is.numeric(n_perm) || length(n_perm) != 1L || is.na(n_perm) ||
      n_perm < 1 || n_perm != round(n_perm) || n_perm > .Machine$integer.max) {
    stop("`n_perm` must be a positive integer.", call. = FALSE)
  }
  as.integer(n_perm)
}

check_stats <- function(t0, t_perm) {
  if (!is.numeric(t0)) stop("`t0` must be a numeric vector.", call. = FALSE)
  t0 <- as.double(t0)
  if (!is.matrix(t_perm) || !is.numeric(t_perm)) {
    stop("`t_perm` must be a numeric matrix: resamples in rows, grid points in columns.",
         call. = FALSE)
  }
  if (ncol(t_perm) != length(t0) || length(t0) < 1L) {
    stop("`t_perm` must have one column per element of `t0`.", call. = FALSE)
  }
  if (nrow(t_perm) < 1L) stop("`t_perm` must have at least one row.", call. = FALSE)
  if (storage.mode(t_perm) != "double") storage.mode(t_perm) <- "double"
  list(t0 = t0, t_perm = t_perm)
}

# Matrices (fd objects are evaluated in prepare.R before reaching here)
check_data <- function(x, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf(paste0("`%s` must be an fd object or a numeric matrix (curves in rows, ",
                        "grid points in columns)."), name), call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(sprintf("`%s` must have at least one curve and one grid point.", name),
         call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop(sprintf("`%s` must not contain NA, NaN or infinite values.", name), call. = FALSE)
  }
  if (storage.mode(x) != "double") storage.mode(x) <- "double"
  x
}

# Numeric mu (fd mu is evaluated in prepare.R)
check_mu <- function(mu, p) {
  if (!is.numeric(mu) || !(length(mu) %in% c(1L, p)) || !all(is.finite(mu))) {
    stop(paste0("`mu` must be a number, a numeric vector with one value per grid point, ",
                "or an fd object (with fd data)."), call. = FALSE)
  }
  mu
}

alt_code <- function(alternative) {
  match(alternative, c("two.sided", "greater", "less")) - 1L
}

# Resampling units: observations (blocks = NULL), q contiguous blocks
# (blocks = q, as cut(seq_len(n), q, labels = FALSE)), or a vector of block ids
# (units in order of first appearance). Returns 0-based ids and their number.
make_units <- function(blocks, n) {
  if (is.null(blocks)) return(list(id = seq_len(n) - 1L, q = n))
  if (length(blocks) == 1L && n > 1L) {
    if (!is.numeric(blocks) || is.na(blocks) || blocks < 1 || blocks > n ||
        blocks != round(blocks)) {
      stop(paste0("A single `blocks` value is the number of contiguous blocks: ",
                  "an integer between 1 and the number of curves."), call. = FALSE)
    }
    id <- if (blocks == 1) rep(1L, n) else cut(seq_len(n), blocks, labels = FALSE)
  } else {
    if (length(blocks) != n) {
      stop(paste0("`blocks` must be NULL, a number of contiguous blocks, or a vector ",
                  "of block ids with one entry per curve (pair)."), call. = FALSE)
    }
    if (anyNA(blocks)) stop("`blocks` must not contain NA.", call. = FALSE)
    id <- blocks
  }
  id <- match(id, unique(id))
  list(id = as.integer(id) - 1L, q = max(id))
}

check_signs <- function(signs, q) {
  if (!is.matrix(signs) || !is.numeric(signs)) {
    stop(paste0("`signs` must be a numeric matrix: one row per resample, one column ",
                "per unit (curve, or block when `blocks` is given)."), call. = FALSE)
  }
  if (ncol(signs) != q) {
    stop(sprintf("`signs` must have %d columns, one per resampling unit.", q), call. = FALSE)
  }
  if (nrow(signs) < 1L || anyNA(signs) || any(signs != 1 & signs != -1)) {
    stop("`signs` must have at least one row and entries -1 or 1.", call. = FALSE)
  }
  signs
}

check_perms <- function(perms, N) {
  if (!is.matrix(perms) || !is.numeric(perms) || ncol(perms) != N || nrow(perms) < 1L) {
    stop(sprintf("`perms` must be a matrix with %d columns: one permutation of 1..%d per row.",
                 N, N), call. = FALSE)
  }
  if (anyNA(perms) || any(perms != round(perms))) {
    stop("`perms` must contain whole numbers only.", call. = FALSE)
  }
  if (any(perms < 1 | perms > N)) {
    stop(sprintf("Each row of `perms` must be a permutation of 1..%d.", N), call. = FALSE)
  }
  storage.mode(perms) <- "integer"
  dimnames(perms) <- NULL
  ok <- apply(perms, 1L, function(r) identical(sort.int(r), seq_len(N)))
  if (!all(ok)) {
    stop(sprintf("Each row of `perms` must be a permutation of 1..%d.", N), call. = FALSE)
  }
  perms
}
