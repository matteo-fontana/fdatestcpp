\name{fdatestcpp-package}
\alias{fdatestcpp-package}
\alias{fdatestcpp}
\docType{package}
\title{Fast Interval-Wise and Threshold-Wise Tests for Functional Data}
\description{
Rcpp port of the parts of the development version of \pkg{fdatest} (Pini,
Stamm and Vantini; \url{https://github.com/permaverse/fdatest}, commit
67bb41c, GPL >= 3) used for domain-selective permutation tests. The CRAN
release 2.1.1 of \pkg{fdatest} provides the earlier interval testing
procedure (\code{ITP*} functions) but not the functions ported here; below,
\pkg{fdatest} refers to that development version.

Three layers:
\itemize{
  \item engine on precomputed statistics: \code{\link{pvalues_pointwise}},
    \code{\link{adjust_iwt}}, \code{\link{adjust_twt}}, \code{\link{adjust_global}};
  \item front end computing observed and resampled pointwise statistics:
    \code{\link{perm_stats1}} (one-sample sign flips) and \code{\link{perm_stats2}}
    (paired swaps, unpaired label permutations), with block-constant signs and
    user-supplied sign or permutation matrices;
  \item drop-in wrappers with \pkg{fdatest}'s argument names, for curves
    given as \code{fd} objects (evaluated on a grid as \pkg{fdatest} does) or
    as matrices:
    \code{\link{iwt1}}, \code{\link{iwt2}}, \code{\link{twt2}}, \code{\link{global2}},
    plus \code{\link{twt1}} and \code{\link{global1}}.
}
The results of the wrappers can be drawn with \code{plot()} and
\code{autoplot()} (\pkg{ggplot2} graphics after the IWTgg scheme, see
\code{\link{autoplot.fdatestcpp_fos}}).
With the same inputs, and the same seed for the wrappers, all p-values are
identical to \pkg{fdatest}'s: the C++ code accumulates in long double in the
same order as R's \code{colMeans()}, \code{rowSums()}, \code{sum()} and
\code{cov()}, and draws the same random numbers. This holds where the
compiler's long double has R's long double precision (as on the usual
x86_64 and arm64 builds; not on R built without long double, where results
agree up to rounding). The deliberate exception is a two-sample \code{mu}
that varies along the grid, which \pkg{fdatest} centres by column by mistake
(see \code{\link{perm_stats2}}).

The functions have \pkg{fdatest}'s names: if both packages are attached, the
one attached last masks the other, so call \code{fdatestcpp::iwt1()} etc.
explicitly when in doubt.
}
\references{
Pini, A. and Vantini, S. (2016). The interval testing procedure: a general
framework for inference in functional data analysis. \emph{Biometrics} 72(3),
835--845. \doi{10.1111/biom.12476}

Pini, A. and Vantini, S. (2017). Interval-wise testing for functional data.
\emph{Journal of Nonparametric Statistics} 29(2), 407--424.
\doi{10.1080/10485252.2017.1306627}

Abramowicz, K., Pini, A., Schelin, L., Sj\enc{ö}{o}stedt de Luna, S., Stamm, A.
and Vantini, S. (2023). Domain selection and familywise error rate for
functional data: a unified framework. \emph{Biometrics} 79(2), 1119--1132.
\doi{10.1111/biom.13669}
}
\keyword{package}
