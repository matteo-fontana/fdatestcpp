# Data preparation: matrices, or fd objects evaluated on a grid exactly as
# fdatest's os_to_coeffs() / ts_to_coeffs() do (same abscissa, same
# t(fda::eval.fd(...)) call), so results stay identical to fdatest's.

is_fd <- function(x) inherits(x, "fd")

# dx is only used (and so only checked) for fd data, as in fdatest; it is
# passed to seq() unchanged
check_dx <- function(dx) {
  if (is.null(dx)) return(NULL)
  if ((!is.numeric(dx) && !is.logical(dx)) || length(dx) != 1L || is.na(dx) ||
      !is.finite(dx) || dx <= 0) {
    stop("`dx` must be NULL or a positive number.", call. = FALSE)
  }
  dx
}

# Grid of an fd object: seq(rangeval[1], rangeval[2], by = dx), dx defaulting
# to 1% of the range (101 points), as in fdatest.
fd_grid <- function(fdobj, dx) {
  dx <- check_dx(dx)
  rangeval <- fdobj$basis$rangeval
  if (is.null(dx)) dx <- (rangeval[2] - rangeval[1]) * 0.01
  # seq() (and so fdatest) fails here too, with a message that does not name dx
  if ((rangeval[2] - rangeval[1]) / dx > .Machine$integer.max) {
    stop("`dx` is too small for the range of the fd object.", call. = FALSE)
  }
  list(rangeval = rangeval, dx = dx, abscissa = seq(rangeval[1], rangeval[2], by = dx))
}

# Curves of an fd object on the abscissa, one per row (fdatest's expression)
fd_values <- function(fdobj, abscissa, name) {
  if (!requireNamespace("fda", quietly = TRUE)) {
    stop("Package 'fda' is needed to evaluate fd objects.", call. = FALSE)
  }
  if (length(dim(fdobj$coefs)) > 2L) {
    stop(sprintf("`%s` is a multivariate fd object (3-d coefficient array): not supported.",
                 name), call. = FALSE)
  }
  x <- t(fda::eval.fd(fdobj = fdobj, evalarg = abscissa))
  if (!all(is.finite(x))) {
    stop(sprintf("`%s` evaluates to NA, NaN or infinite values on the grid.", name),
         call. = FALSE)
  }
  x
}

# mu: a number, a vector with one value per grid point, or an fd object on the
# same range as fd data (evaluated on the data's grid; fdatest returns the
# 1 x p matrix of its values as `mu`).
prep_mu <- function(mu, p, grid) {
  if (is_fd(mu)) {
    if (is.null(grid)) {
      stop(paste0("`mu` is an fd object but the data are matrices: evaluate `mu` on ",
                  "the data's grid and pass the vector of values."), call. = FALSE)
    }
    rangeval_mu <- mu$basis$rangeval
    if (sum(rangeval_mu == grid$rangeval) != 2) {
      stop("The range of values of `mu` must be the same as the range of values of the data.",
           call. = FALSE)
    }
    mu_eval <- fd_values(mu, seq(rangeval_mu[1], rangeval_mu[2], by = grid$dx), "mu")
    if (nrow(mu_eval) != 1L) stop("`mu` must be a single function.", call. = FALSE)
    return(mu_eval)
  }
  check_mu(mu, p)
}

# One sample: list(data = curves as fdatest reports them, x = numeric matrix,
# mu = mu as fdatest reports it, abscissa = grid of fd input or NULL)
prep_one <- function(data, mu, dx) {
  if (is_fd(data)) {
    grid <- fd_grid(data, dx)
    x <- fd_values(data, grid$abscissa, "data")
  } else {
    grid <- NULL
    x <- check_data(data, "data")
  }
  list(data = if (is.matrix(data)) data else x, x = x,
       mu = prep_mu(mu, ncol(x), grid), abscissa = grid$abscissa)
}

# Two samples: both fd objects on the same range, or both matrices
prep_two <- function(data1, data2, mu, dx) {
  if (is_fd(data1) && is_fd(data2)) {
    if (sum(data1$basis$rangeval == data2$basis$rangeval) != 2) {
      stop(paste0("The range of values of `data1` must be the same as the range of ",
                  "values of `data2`."), call. = FALSE)
    }
    grid <- fd_grid(data1, dx)
    x1 <- fd_values(data1, grid$abscissa, "data1")
    x2 <- fd_values(data2, grid$abscissa, "data2")
  } else if (!is_fd(data1) && !is_fd(data2)) {
    grid <- NULL
    x1 <- check_data(data1, "data1")
    x2 <- check_data(data2, "data2")
  } else {
    stop(paste0("Both `data1` and `data2` must be either functional data objects of ",
                "class fd or matrices."), call. = FALSE)
  }
  if (ncol(x2) != ncol(x1)) {
    stop("`data1` and `data2` must have the same number of columns (grid points).",
         call. = FALSE)
  }
  list(data = rbind(if (is.matrix(data1)) data1 else x1, if (is.matrix(data2)) data2 else x2),
       x1 = x1, x2 = x2, mu = prep_mu(mu, ncol(x1), grid), abscissa = grid$abscissa)
}
