\name{iwt1}
\alias{iwt1}
\alias{twt1}
\alias{global1}
\title{One-sample interval-wise, threshold-wise and global tests}
\description{
Drop-in version of \code{iwt1()} from the development version of
\pkg{fdatest} (see \code{\link{fdatestcpp-package}}; here \pkg{fdatest} means
that version, as the CRAN release 2.1.1 has no \code{iwt1()}): tests the
center of symmetry \code{mu} of a functional sample by sign flips, with
interval-wise adjusted p-values. \code{twt1()} and \code{global1()} (not in
\pkg{fdatest}) give threshold-wise adjusted and global p-values on the same
resamples.
}
\usage{
iwt1(data, mu = 0, n_perm = 1000L, dx = NULL, recycle = TRUE, verbose = FALSE,
     aggregation_strategy = c("integral", "max"),
     alternative = c("two.sided", "less", "greater"), standardize = FALSE,
     blocks = NULL, signs = NULL, plus_one = FALSE, return_perm = FALSE)

twt1(data, mu = 0, n_perm = 1000L, dx = NULL, verbose = FALSE,
     aggregation_strategy = c("integral", "max"),
     alternative = c("two.sided", "less", "greater"), standardize = FALSE,
     blocks = NULL, signs = NULL, plus_one = FALSE, return_perm = FALSE)

global1(data, mu = 0, n_perm = 1000L, dx = NULL, verbose = FALSE,
        aggregation_strategy = c("integral", "max"),
        alternative = c("two.sided", "less", "greater"), standardize = FALSE,
        blocks = NULL, signs = NULL, plus_one = FALSE, return_perm = FALSE)
}
\arguments{
\item{data, mu, n_perm, dx, alternative, standardize, blocks, signs}{see
  \code{\link{perm_stats1}}: \code{data} is an \code{fd} object or a matrix,
  \code{dx} the step of the grid on which \code{fd} objects are evaluated.}
\item{verbose}{accepted for compatibility with \pkg{fdatest} and ignored.}
\item{recycle, aggregation_strategy, plus_one}{see \code{\link{adjust_iwt}}.
  \code{recycle} defaults to \code{TRUE} as in \pkg{fdatest}: set it to
  \code{FALSE} unless the domain is periodic.}
\item{return_perm}{also return the \code{B x p} matrix of resampled statistics.}
}
\value{
An object of class \code{c("fdatestcpp_fos", "fos")} with \pkg{fdatest}'s fields \code{data},
\code{mu}, \code{unadjusted_pvalues}, \code{adjusted_pvalues},
\code{correction_method} and, for the IWT, \code{pvalue_matrix}; plus
\code{global_pvalue}, \code{t0}, if requested \code{t_perm}, and for
\code{fd} data \code{abscissa}, the grid on which the p-values are given.
As in \pkg{fdatest}, \code{data} holds the curves evaluated on the grid and
an \code{fd} \code{mu} is returned as the \code{1 x p} matrix of its values.
\code{plot()} and \code{autoplot()} draw the result with \pkg{ggplot2}
(see \code{\link{autoplot.fdatestcpp_fos}}).
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
\code{\link{iwt2}} for two-sample tests, \code{\link{perm_stats1}} and
\code{\link{adjust_iwt}} for the front end and the engine used here.
}
\examples{
set.seed(1)
d <- matrix(rnorm(40 * 25), 40, 25)
d[, 18:25] <- d[, 18:25] + 0.7
r <- iwt1(d, n_perm = 500, recycle = FALSE)
which(r$adjusted_pvalues <= 0.05)
twt1(d, n_perm = 500, blocks = 10)$adjusted_pvalues

# fd input, evaluated on 101 points of its range
tt <- seq(0, 1, length.out = 25)
f <- fda::smooth.basis(tt, t(d), fda::create.bspline.basis(c(0, 1), 10))$fd
r <- iwt1(f, n_perm = 500, recycle = FALSE)
r$abscissa[r$adjusted_pvalues <= 0.05]
global1(f, n_perm = 500)$global_pvalue
}
