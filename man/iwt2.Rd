\name{iwt2}
\alias{iwt2}
\alias{twt2}
\alias{global2}
\title{Two-sample interval-wise, threshold-wise and global tests}
\description{
Drop-in versions of \code{iwt2()}, \code{twt2()} and \code{global2()} from
the development version of \pkg{fdatest} (see
\code{\link{fdatestcpp-package}}; here \pkg{fdatest} means that version, as
the CRAN release 2.1.1 lacks these functions): test the mean difference
\code{mu} between two functional samples, paired (swaps within pairs) or not
(label permutations), with interval-wise, threshold-wise or global adjusted
p-values.
}
\usage{
iwt2(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
     alternative = c("two.sided", "less", "greater"), standardize = FALSE,
     verbose = FALSE, aggregation_strategy = c("integral", "max"),
     recycle = TRUE, blocks = NULL, signs = NULL, perms = NULL,
     plus_one = FALSE, return_perm = FALSE)

twt2(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
     alternative = c("two.sided", "less", "greater"), standardize = FALSE,
     verbose = FALSE, aggregation_strategy = c("integral", "max"),
     blocks = NULL, signs = NULL, perms = NULL,
     plus_one = FALSE, return_perm = FALSE)

global2(data1, data2, mu = 0, dx = NULL, n_perm = 1000L, paired = FALSE,
        alternative = c("two.sided", "less", "greater"), standardize = FALSE,
        verbose = FALSE, aggregation_strategy = c("integral", "max"),
        blocks = NULL, signs = NULL, perms = NULL,
        plus_one = FALSE, return_perm = FALSE)
}
\arguments{
\item{data1, data2, mu, dx, n_perm, paired, alternative, standardize, blocks,
  signs, perms}{see \code{\link{perm_stats2}}: \code{data1}, \code{data2} are
  both \code{fd} objects on the same range or both matrices; \code{dx} is the
  step of the grid on which \code{fd} objects are evaluated.}
\item{verbose}{accepted for compatibility with \pkg{fdatest} and ignored.}
\item{recycle, aggregation_strategy, plus_one}{see \code{\link{adjust_iwt}}.
  \code{recycle} defaults to \code{TRUE} as in \pkg{fdatest}: set it to
  \code{FALSE} unless the domain is periodic.}
\item{return_perm}{also return the \code{B x p} matrix of resampled statistics.}
}
\value{
An object of class \code{c("fdatestcpp_fts", "fts")} with \pkg{fdatest}'s fields \code{data},
\code{group_labels}, \code{mu}, \code{unadjusted_pvalues},
\code{adjusted_pvalues}, \code{correction_method} and, for the IWT,
\code{pvalue_matrix}; plus \code{global_pvalue}, \code{t0}, if requested
\code{t_perm}, and for \code{fd} data \code{abscissa}, the grid on which the
p-values are given. As in \pkg{fdatest}, \code{data} holds the curves
evaluated on the grid and an \code{fd} \code{mu} is returned as the
\code{1 x p} matrix of its values (\code{global2()}'s \code{global_pvalue} is
\pkg{fdatest}'s own field). \code{plot()} and \code{autoplot()} draw the
result with \pkg{ggplot2} (see \code{\link{autoplot.fdatestcpp_fts}}).
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
\code{\link{iwt1}} for one-sample tests, \code{\link{perm_stats2}} and
\code{\link{adjust_iwt}} for the front end and the engine used here.
}
\examples{
set.seed(1)
l1 <- matrix(rnorm(48 * 20), 48, 20)
l2 <- l1 + matrix(rnorm(48 * 20, sd = 0.5), 48, 20)
l2[, 15:20] <- l2[, 15:20] - 0.5
r <- iwt2(l1, l2, paired = TRUE, n_perm = 500, recycle = FALSE)
round(r$adjusted_pvalues, 3)
# swaps constant within 12 blocks of consecutive pairs
twt2(l1, l2, paired = TRUE, n_perm = 500, blocks = 12)$adjusted_pvalues
global2(l1, l2, paired = TRUE, n_perm = 500)$global_pvalue
}
