# Front end: observed and resampled statistics.
# With the same seed, t0 and t_perm must be identical() to fdatest's
# os_permtest() / ts_permtest(); blocks, signs and perms are checked against
# direct R computations.

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
  all(vapply(c("os_permtest", "ts_permtest"), exists, NA,
             envir = asNamespace("fdatest"), inherits = FALSE))

n_checks <- 0L
fails <- character(0)
expect <- function(ok, label) {
  n_checks <<- n_checks + 1L
  if (!isTRUE(ok)) fails <<- c(fails, label)
}

set.seed(1)
n <- 25L
p <- 13L
B <- 200L
tau <- seq_len(p) / p
R <- chol(exp(-abs(outer(tau, tau, "-")) / 0.3))
curves <- function(n, shift = 0) matrix(rnorm(n * p), n, p) %*% R + shift
x1 <- curves(n, 0.3)
x2 <- curves(n) + 0.5 * x1      # paired with x1
x3 <- curves(n + 7L)            # unpaired, different size

if (has_fdatest) {
  fd <- asNamespace("fdatest")
  for (mu in list(0, 0.2, seq(-0.1, 0.1, length.out = p))) {
    coeff <- x1 - matrix(mu, nrow = n, ncol = p, byrow = TRUE)
    set.seed(5); a <- fd$os_permtest(coeff, B)
    set.seed(5); b <- perm_stats1(x1, mu = mu, n_perm = B)
    lab <- sprintf("one-sample mu length %d", length(mu))
    expect(identical(b$t0, a$t0), paste(lab, "t0"))
    expect(identical(b$t_perm, a$t_coeff), paste(lab, "t_perm"))
  }
  for (paired in c(TRUE, FALSE)) {
    for (alt in c("two.sided", "greater", "less")) {
      for (std in c(FALSE, TRUE)) {
        d2 <- if (paired) x2 else x3
        coeff <- rbind(x1 - 0.1, d2)
        set.seed(7); a <- fd$ts_permtest(coeff, n, B, alt, paired, std)
        set.seed(7); b <- perm_stats2(x1, d2, mu = 0.1, n_perm = B, paired = paired,
                                      alternative = alt, standardize = std)
        lab <- sprintf("two-sample paired=%s %s standardize=%s", paired, alt, std)
        expect(identical(b$t0, a$t0), paste(lab, "t0"))
        expect(identical(b$t_perm, a$t_coeff), paste(lab, "t_perm"))
        set.seed(7); u1 <- runif(1)
        set.seed(7); invisible(perm_stats2(x1, d2, mu = 0.1, n_perm = B, paired = paired,
                                           alternative = alt, standardize = std))
        u2 <- runif(1)
        set.seed(7); invisible(fd$ts_permtest(coeff, n, B, alt, paired, std))
        u3 <- runif(1)
        expect(identical(u2, u3) && !identical(u1, u2), paste(lab, "RNG stream"))
      }
    }
  }
}

# User-supplied signs / perms reproduce the internal draws
set.seed(9); a <- perm_stats1(x1, n_perm = B)
set.seed(9); S <- matrix(rbinom(B * n, 1, 0.5), B, n, byrow = TRUE) * 2 - 1
expect(identical(perm_stats1(x1, signs = S)$t_perm, a$t_perm), "signs, one-sample")
set.seed(9); a <- perm_stats2(x1, x2, n_perm = B, paired = TRUE, alternative = "less")
set.seed(9); W <- matrix(rbinom(B * n, 1, 0.5), B, n, byrow = TRUE)
b <- perm_stats2(x1, x2, paired = TRUE, alternative = "less", signs = 1 - 2 * W)
expect(identical(b$t_perm, a$t_perm), "signs, paired (-1 = swap)")
set.seed(9); a <- perm_stats2(x1, x3, n_perm = B, standardize = TRUE)
set.seed(9); P <- t(replicate(B, sample(2L * n + 7L)))
expect(identical(perm_stats2(x1, x3, perms = P, standardize = TRUE)$t_perm, a$t_perm),
       "perms, unpaired")

# Blocks
set.seed(4); a <- perm_stats1(x1, n_perm = B)
set.seed(4); b <- perm_stats1(x1, n_perm = B, blocks = seq_len(n))
expect(identical(a$t_perm, b$t_perm) && identical(a$t0, b$t0), "singleton blocks = no blocks")
q <- 6L
id <- cut(seq_len(n), q, labels = FALSE)
set.seed(5); Sb <- matrix(sample(c(-1, 1), B * q, TRUE), B, q)
manual <- function(f) t(apply(Sb, 1, function(s) f(x1 * s[id])))
b <- perm_stats1(x1, signs = Sb, blocks = q)
expect(isTRUE(all.equal(b$t_perm, manual(function(y) colMeans(y)^2), tolerance = 1e-13)),
       "block flips vs R")
expect(isTRUE(all.equal(b$t0, colMeans(x1)^2, tolerance = 1e-14)), "block t0 vs R")
expect(identical(perm_stats1(x1, signs = Sb, blocks = id * 10)$t_perm, b$t_perm),
       "block id vector = block count")
expect(b$n_units == q, "number of units")
b <- perm_stats1(x1, signs = Sb, blocks = q, standardize = TRUE)
# (identical() on x86_64 Windows/Linux; a tolerance keeps the test portable to
# platforms where compilers may round R's and the package's code differently)
expect(isTRUE(all.equal(b$t_perm, manual(function(y) colMeans(y)^2 / apply(y, 2, var)),
                        tolerance = 1e-13)),
       "standardized block flips = R's colMeans(y)^2 / apply(y, 2, var)")
b <- perm_stats1(x1, signs = Sb, blocks = q, alternative = "greater")
expect(isTRUE(all.equal(b$t_perm, manual(function(y) pmax(colMeans(y), 0)^2),
                        tolerance = 1e-13)), "greater block flips vs R")
b <- perm_stats1(x1, signs = Sb, blocks = q, alternative = "less")
expect(isTRUE(all.equal(b$t_perm, manual(function(y) pmin(colMeans(y), 0)^2),
                        tolerance = 1e-13)), "less block flips vs R")
bp <- perm_stats2(x1, x2, paired = TRUE, signs = Sb, blocks = q)
bo <- perm_stats1(x1 - x2, signs = Sb, blocks = q)
expect(isTRUE(all.equal(bp$t_perm, bo$t_perm, tolerance = 1e-12)),
       "paired block swaps = one-sample block flips of the difference")
# paired standardized with blocks: swaps constant within blocks, fdatest's statistic
bs <- perm_stats2(x1, x2, paired = TRUE, signs = Sb, blocks = q, standardize = TRUE)
man <- t(apply(Sb, 1, function(s) {
  sw <- s[id] < 0
  g1 <- x1; g1[sw, ] <- x2[sw, ]
  g2 <- x2; g2[sw, ] <- x1[sw, ]
  (colMeans(g1) - colMeans(g2))^2 / ((apply(g1, 2, var) + apply(g2, 2, var)) / 2)
}))
expect(isTRUE(all.equal(bs$t_perm, man, tolerance = 1e-12)), "paired standardized blocks vs R")

# argument checks
err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
expect(err(perm_stats1(x1, signs = S[, -1])), "signs width")
expect(err(perm_stats1(x1, signs = S * 2)), "signs values")
expect(err(perm_stats2(x1, x3, paired = TRUE)), "paired sizes")
expect(err(perm_stats2(x1, x3, blocks = 3)), "blocks unpaired")
expect(err(perm_stats1(x1, blocks = rep(1:2, length.out = n - 1))), "blocks length")
expect(err(perm_stats1(rbind(x1, NA))), "non-finite data")

cat(sprintf("test-resample: %d checks, %d failed%s\n", n_checks, length(fails),
            if (has_fdatest) "" else " (fdatest comparisons skipped)"))
if (length(fails)) stop("Failed:\n", paste(fails, collapse = "\n"))
