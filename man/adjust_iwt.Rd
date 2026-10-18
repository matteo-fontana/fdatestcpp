\name{adjust_iwt}
\alias{adjust_iwt}
\alias{adjust_twt}
\alias{adjust_global}
\alias{pvalues_pointwise}
\title{P-values from observed and resampled pointwise statistics}
\description{
The engine: from observed pointwise statistics \code{t0} and resampled ones
\code{t_perm}, compute unadjusted p-values and interval-wise (IWT),
threshold-wise (TWT) or global adjusted p-values exactly as the internal
functions \code{p_adjust_iwt()}, \code{p_adjust_twt()} and
\code{p_adjust_global()} of the development version of \pkg{fdatest} (see
\code{\link{fdatestcpp-package}}; here and below \pkg{fdatest} means that
version). Any resampling scheme can feed it (block sign flips, sieve or wild
bootstraps, ...).
}
\usage{
pvalues_pointwise(t0, t_perm, plus_one = FALSE)

adjust_iwt(t0, t_perm, recycle, aggregation_strategy = c("integral", "max"),
           plus_one = FALSE, pvalue_matrix = TRUE)

adjust_twt(t0, t_perm, aggregation_strategy = c("integral", "max"),
           plus_one = FALSE)

adjust_global(t0, t_perm, aggregation_strategy = c("integral", "max"),
              plus_one = FALSE)
}
\arguments{
\item{t0}{numeric vector of length \code{p}: observed pointwise statistics
  (large values speak against the null).}
\item{t_perm}{numeric \code{B x p} matrix: resampled pointwise statistics, one
  resample per row.}
\item{recycle}{\code{TRUE} for a periodic domain (intervals wrap around),
  \code{FALSE} otherwise. Required: \pkg{fdatest}'s \code{iwt1()}/\code{iwt2()}
  default to \code{TRUE}.}
\item{aggregation_strategy}{how pointwise statistics are combined over an
  interval or set: \code{"integral"} (sum) or \code{"max"}.}
\item{plus_one}{\code{FALSE} (default, as \pkg{fdatest}): p-value =
  \code{mean(T* >= T0)}. \code{TRUE}: \code{(1 + #(T* >= T0)) / (B + 1)}.}
\item{pvalue_matrix}{logical: if \code{TRUE} (default), \code{adjust_iwt()}
  also returns the \code{p x p} matrix of interval p-values.}
}
\details{
IWT: every interval's statistic is accumulated from its start, one point at a
time (a long double running sum rounded to double, or a running maximum), so
all \eqn{p(p+1)/2} intervals (\eqn{p^2} when recycling) cost \eqn{O(B p^2)} and
reproduce R's \code{sum()}/\code{rowSums()} rounding. Adjusted p-values are the
maximum interval p-value over intervals containing each point, computed by an
\eqn{O(p^2)} recursion.

TWT: grid points are sorted by unadjusted p-value, so the family of sub- and
super-level sets of the p-value function becomes prefixes and suffixes of that
order and all set statistics of a resample come out of one pass, \eqn{O(B p)}.
\pkg{fdatest} sums each set in grid order; when a sorted-order sum is within
its rounding-error bound of the observed value the set is re-summed in grid
order, so the comparisons, and the p-values, are identical.

\code{pvalue_matrix} follows \pkg{fdatest}: row \code{p} holds the unadjusted
p-values; row \code{i} holds intervals of length \code{p - i + 1}, entry
\code{[i, j]} the interval starting at \code{j}; with \code{recycle = FALSE}
row \code{i} has entries \code{j = 1..i} and \code{NA} elsewhere.

Non-finite statistics are handled as \pkg{fdatest} does: \code{NaN} is dropped
inside IWT interval and TWT set statistics (\code{na.rm = TRUE}) and makes the
unadjusted p-value of its grid point \code{NA}; \code{adjust_global()} sums
without \code{na.rm}, like \code{fdatest:::p_adjust_global()}, so any
\code{NaN} statistic makes its p-values \code{NA}. Sums that overflow are also
treated as \pkg{fdatest} treats them: observed sums use \code{sum()}'s rule
(beyond the double range: \code{Inf}), resampled sums \code{rowSums()}'s
(rounded to nearest, so totals below \code{DBL_MAX + 2^970} give
\code{DBL_MAX}); TWT sets whose sums approach overflow are decided by the
exact grid-order sums.

As in \pkg{fdatest}, the unadjusted p-values (and the TWT adjusted p-values)
are named after the column names of \code{t_perm}, if any.
}
\value{
\code{pvalues_pointwise()} returns the numeric vector of length \code{p} of
unadjusted p-values. The other functions return a list with
\item{unadjusted_pvalues}{numeric vector of length \code{p}: the unadjusted
  p-values.}
\item{adjusted_pvalues}{numeric vector of length \code{p}: the adjusted
  p-value function (constant for \code{adjust_global()}).}
\item{pvalue_matrix}{\code{adjust_iwt()} with \code{pvalue_matrix = TRUE}
  only: \code{p x p} numeric matrix of interval p-values in \pkg{fdatest}'s
  layout (see Details).}
\item{global_pvalue}{a number, the p-value of the whole domain:
  \code{pvalue_matrix[1, 1]} for the IWT; the full level set for the TWT; the
  global test for \code{adjust_global()}.}
}
\references{
Pini, A. and Vantini, S. (2017). Interval-wise testing for functional data.
\emph{Journal of Nonparametric Statistics} 29(2), 407--424.
\doi{10.1080/10485252.2017.1306627}

Abramowicz, K., Pini, A., Schelin, L., Sj\enc{ö}{o}stedt de Luna, S., Stamm, A.
and Vantini, S. (2023). Domain selection and familywise error rate for
functional data: a unified framework. \emph{Biometrics} 79(2), 1119--1132.
\doi{10.1111/biom.13669}
}
\seealso{
\code{\link{perm_stats1}} and \code{\link{perm_stats2}} to compute
\code{t0} and \code{t_perm}; \code{\link{iwt1}} and \code{\link{iwt2}} for the
complete tests.
}
\examples{
set.seed(1)
x <- matrix(rnorm(30 * 20), 30, 20) + rep(c(0, 0.6), each = 300)
t0 <- colMeans(x)^2
t_perm <- t(replicate(500, colMeans(x * sample(c(-1, 1), 30, TRUE))^2))
round(pvalues_pointwise(t0, t_perm), 3)
r <- adjust_iwt(t0, t_perm, recycle = FALSE)
round(r$adjusted_pvalues, 3)
r$global_pvalue
adjust_twt(t0, t_perm)$adjusted_pvalues
adjust_global(t0, t_perm)$global_pvalue
}
