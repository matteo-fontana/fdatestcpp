\name{perm_stats1}
\alias{perm_stats1}
\alias{perm_stats2}
\title{Observed and resampled pointwise statistics}
\description{
The front end: pointwise test statistics of the data and of \code{n_perm}
resamples, ready for \code{\link{adjust_iwt}}, \code{\link{adjust_twt}} or
\code{\link{adjust_global}}.

\code{perm_stats1()}: one-sample test of the center of symmetry \code{mu},
resampling by sign flips of the centered curves (as \code{iwt1()} of
\pkg{fdatest}'s development version; here \pkg{fdatest} means that version,
see \code{\link{fdatestcpp-package}}). \code{perm_stats2()}: two-sample test
of the mean difference \code{mu}, by swapping the curves within pairs
(\code{paired = TRUE}, the same as sign flips of the differences) or by
permuting the labels.
}
\usage{
perm_stats1(data, mu = 0, n_perm = 1000L, dx = NULL,
            alternative = c("two.sided", "less", "greater"),
            standardize = FALSE, blocks = NULL, signs = NULL)

perm_stats2(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
            alternative = c("two.sided", "less", "greater"),
            standardize = FALSE, blocks = NULL, signs = NULL, perms = NULL)
}
\arguments{
\item{data, data1, data2}{functional data objects of class \code{fd}
  (\pkg{fda}), or numeric matrices with one curve per row evaluated on a
  common grid of \code{p} points (columns). Two-sample data must be both
  \code{fd} objects on the same range or both matrices.}
\item{mu}{null value: a number, a vector with one value per grid point, or
  (with \code{fd} data) an \code{fd} object with one function on the same
  range, evaluated on the data's grid.}
\item{dx}{step of the grid on which \code{fd} objects are evaluated:
  \code{seq(range[1], range[2], by = dx)}, as \pkg{fdatest}; \code{NULL} (the
  default) uses 1\% of the range, i.e. 101 points. Ignored (and not checked)
  for matrices, as in \pkg{fdatest}.}
\item{n_perm}{number of resamples \code{B} (ignored when \code{signs} or
  \code{perms} is given).}
\item{paired}{paired swaps (\code{TRUE}) or label permutations.}
\item{alternative}{\code{"two.sided"}: squared mean (difference);
  \code{"greater"} / \code{"less"}: the same, set to 0 when its sign is not the
  alternative's. As \pkg{fdatest}'s two-sample tests; one-sample tests apply
  it to the mean of \code{data - mu}.}
\item{standardize}{divide by a variance recomputed on every resample.
  Two-sample: \pkg{fdatest}'s pooled variance of the two (swapped or permuted)
  samples. One-sample (not in \pkg{fdatest}): \code{var()} of the
  sign-flipped curves \code{y}, so that the statistic is exactly R's
  \code{colMeans(y)^2 / apply(y, 2, var)}: \code{+Inf} where all curves are
  equal to a nonzero value, \code{NaN} for 0/0 (all curves zero at a grid
  point, or the wrong sign under a one-sided alternative). \code{NaN}
  statistics give \code{NA} unadjusted p-values, are dropped by the IWT and
  TWT and make the Global test's p-values \code{NA}, as in \pkg{fdatest}.}
\item{blocks}{\code{NULL} (every curve or pair is its own unit), a number
  \code{q} of contiguous blocks (\code{cut(seq_len(n), q, labels = FALSE)}), or
  a vector of block ids, one per curve or pair, whose units are ordered by
  first appearance. Signs (swaps) are constant within a block. Sign flips and
  paired swaps only.}
\item{signs}{optional \code{B x q} matrix of \code{-1}/\code{1}, one column
  per unit (curve, pair or block; for a vector of block ids, column \code{g}
  is the \code{g}-th id to appear, \code{unique(blocks)[g]}). One-sample:
  multiplies the unit's centered curves. Paired: \code{-1} swaps the unit's
  curves.}
\item{perms}{unpaired only: optional \code{B x (n1 + n2)} matrix whose rows are
  permutations of \code{1:(n1 + n2)}, applied to \code{rbind(data1, data2)};
  the first \code{n1} rows of the permuted data form sample 1.}
}
\details{
Random signs are drawn as \code{rbinom(q, 1, 0.5)} per resample and
permutations as \code{sample(n1 + n2)}, in \pkg{fdatest}'s order, so with
\code{blocks = NULL} the same seed gives the same resamples and identical
statistics as \pkg{fdatest}. One-sample draws of 1 keep the sign; paired draws
of 1 swap the pair. With \code{signs}, row \code{b} defines resample \code{b}.

With blocks, a unit's curves are summed once and the statistics of each
resample cost \eqn{O(q p)} instead of \eqn{O(n p)}, except with
\code{standardize = TRUE}, which is computed curve by curve. Built from block
sums, \code{t0} can then differ from the unblocked statistic in the last bits;
the identity resample reproduces it exactly.

\code{fd} objects are evaluated exactly as in \pkg{fdatest}
(\code{t(fda::eval.fd(...))} on the grid above), so their statistics are
identical to \pkg{fdatest}'s.

Two-sample \code{mu} is subtracted from column \code{j} of \code{data1};
\pkg{fdatest} subtracts \code{matrix(mu, n1, p)}, which fills by column. The
two agree only when the values of \code{mu} on the grid are exactly equal: a
scalar, a vector of equal values, or an \code{fd} constant that evaluates to
exactly equal values (for example on a constant basis, or the constant term
of a Fourier basis). A constant
function on a B-spline basis evaluates to values that differ in the last
bits, so its statistics differ from \pkg{fdatest}'s at that level.
}
\value{
A list with \code{t0} (length \code{p}), \code{t_perm} (\code{B x p}),
\code{n_units} (the number of resampling units), \code{mu} (as used: a
number, a vector, or the \code{1 x p} matrix of an \code{fd} \code{mu}'s
values) and, for \code{fd} data, \code{abscissa} (the grid). As in
\pkg{fdatest}, \code{t0} and the columns of \code{t_perm} are named after the
data's column names, if any.
}
\examples{
set.seed(1)
d <- matrix(rnorm(48 * 15), 48, 15)
s <- perm_stats1(d, n_perm = 200, blocks = 12)   # signs constant within 12 blocks
dim(s$t_perm)
adjust_iwt(s$t0, s$t_perm, recycle = FALSE)$global_pvalue

# user-supplied signs, e.g. from another scheme
S <- matrix(sample(c(-1, 1), 200 * 12, TRUE), 200, 12)
s2 <- perm_stats1(d, signs = S, blocks = 12)
adjust_twt(s2$t0, s2$t_perm)$global_pvalue

# fd objects are evaluated on seq(range[1], range[2], by = dx)
tt <- seq(0, 1, length.out = 15)
f <- fda::smooth.basis(tt, t(d), fda::create.bspline.basis(c(0, 1), 8))$fd
s3 <- perm_stats1(f, n_perm = 200, dx = 0.05)
length(s3$abscissa)

# two samples: paired swaps, and label permutations given by the user
s4 <- perm_stats2(d[1:24, ], d[25:48, ] + 0.5, n_perm = 200, paired = TRUE)
adjust_iwt(s4$t0, s4$t_perm, recycle = FALSE)$global_pvalue
P <- t(replicate(200, sample(48)))
s5 <- perm_stats2(d[1:24, ], d[25:48, ], perms = P)
dim(s5$t_perm)
}
