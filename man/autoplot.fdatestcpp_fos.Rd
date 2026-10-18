\name{autoplot.fdatestcpp_fos}
\alias{autoplot.fdatestcpp_fos}
\alias{autoplot.fdatestcpp_fts}
\alias{plot.fdatestcpp_fos}
\alias{plot.fdatestcpp_fts}
\title{Plot the Result of an Interval-Wise, Threshold-Wise or Global Test}
\description{
\pkg{ggplot2} graphics for the results of \code{\link{iwt1}},
\code{\link{twt1}}, \code{\link{global1}}, \code{\link{iwt2}},
\code{\link{twt2}} and \code{\link{global2}}, after the plotting scheme of
the IWTgg code by Matteo Fontana: the curves, and the adjusted and
unadjusted p-value functions, in two panels over a common abscissa, with the
part of the domain where the adjusted p-value is below \code{alpha1} shaded
in light grey, and below \code{alpha2} in darker grey.
}
\usage{
\method{autoplot}{fdatestcpp_fos}(object, panels = c("data", "pvalues"),
  alpha1 = 0.05, alpha2 = 0.01, xrange = NULL, unadjusted = TRUE,
  steps = NULL, means = TRUE, title = NULL, xlab = "Domain",
  ylabel = "Functional data", linewidth = 0.5, ...)

\method{autoplot}{fdatestcpp_fts}(object, panels = c("data", "pvalues"),
  alpha1 = 0.05, alpha2 = 0.01, xrange = NULL, unadjusted = TRUE,
  steps = NULL, means = TRUE, title = NULL, xlab = "Domain",
  ylabel = "Functional data", linewidth = 0.5, ...)

\method{plot}{fdatestcpp_fos}(x, ...)

\method{plot}{fdatestcpp_fts}(x, ...)
}
\arguments{
\item{object, x}{the result of \code{iwt1()}, \code{twt1()}, \code{global1()}
  (class \code{"fdatestcpp_fos"}) or of \code{iwt2()}, \code{twt2()},
  \code{global2()} (class \code{"fdatestcpp_fts"}).}
\item{panels}{which panels to draw, from top to bottom: \code{"data"} (the
  curves) and \code{"pvalues"} (the p-value functions).}
\item{alpha1, alpha2}{the two significance levels shaded in grey, in
  \code{(0, 1]}; the larger is shaded lighter. As in \pkg{fdatest} and
  IWTgg, a point is shaded when its adjusted p-value is strictly below the
  level, and the two values are swapped if \code{alpha1 < alpha2}.}
\item{xrange}{\code{NULL}, or the range of the abscissa as in
  \pkg{fdatest}'s plots: the p-values are then drawn at \code{p} equally
  spaced points between \code{xrange[1]} and \code{xrange[2]}. With
  \code{NULL}, the grid of \code{fd} data (\code{object$abscissa}) is used,
  and \code{[0, 1]} for matrix data. For \code{fd} data the default is
  exact: an \code{xrange} that moves the points off the evaluation grid (for
  example the basis range when \code{dx} does not divide it) gives a
  warning. To relabel the axis of \code{fd} data, transform
  \code{object$abscissa} instead.}
\item{unadjusted}{also draw the unadjusted p-values (dashed).}
\item{steps}{draw the p-value functions as steps, constant over the cell of
  each grid point (the cells that are shaded), instead of lines joining the
  grid points. \code{NULL} (default) uses steps on coarse grids, of at most
  50 points, where a line between two grid points would cross a level away
  from the edge of the shading; \code{TRUE} or \code{FALSE} forces either.}
\item{means}{draw the mean curve of each sample (thick lines) and the
  null hypothesis (dashed): \code{mu} in one-sample tests; in two-sample
  tests with \code{mu} not 0, the mean of the second sample plus \code{mu},
  where the mean of the first sample lies under the null hypothesis (drawn
  in the first sample's colour and explained in the plot caption).}
\item{title}{plot title.}
\item{xlab}{title of the x axis.}
\item{ylabel}{label of the data panel, a single string (\pkg{fdatest}'s
  argument name).}
\item{linewidth}{width of the curves; means and p-value functions are drawn
  wider.}
\item{...}{not used by \code{autoplot()}, which warns about any argument
  given here (such as base-graphics \code{main} or \code{lwd}). \code{plot()}
  passes its arguments to \code{autoplot()}.}
}
\details{
The panels are facets of a single ggplot, so the result can be modified with
the usual \pkg{ggplot2} functions and saved with \code{ggplot2::ggsave()}.
For example, \code{+ ggplot2::scale_colour_manual(values = c("1" = "red",
"2" = "darkgreen"), name = "Sample")} recolours the two samples, their means
and the dashed reference line. The facet strips on the left serve as y-axis titles and the
legends are stacked below the plot; a complete theme such as
\code{ggplot2::theme_bw()} resets both, so follow it with
\code{ggplot2::theme(strip.placement = "outside", legend.position = "bottom",
legend.box = "vertical", legend.location = "plot")} to keep the layout.

The two samples are drawn in the colour-blind-safe orange and blue of the
Okabe-Ito palette, interleaved so that neither covers the other, with
transparency that grows with the number of curves (the lighter orange
curves more opaque, so that both samples stand out equally); the means are
drawn over a white halo. The legends come in a fixed order: shading,
curves, p-values. The grid is drawn over the shading.
Each grid point is shaded over the cell between the midpoints to its
neighbours, so that a single significant point is visible and adjacent
points form one band. With a single grid point, curves and p-values are drawn
as points.

For \code{global1()} and \code{global2()} the adjusted p-value is the global
p-value at every point, so the whole domain is shaded when the global test
rejects.

The objects keep \pkg{fdatest}'s classes \code{"fos"} and \code{"fts"}
after their own, so \pkg{fdatest}'s \code{IWTimage()} also works on IWT
results; \pkg{fdatest}'s \code{plot()} and \code{autoplot()} methods can be
called explicitly as \code{fdatest:::plot.fos(r)}.
}
\value{
\code{autoplot()} returns a ggplot object. \code{plot()} draws it and
returns it invisibly.
}
\seealso{\code{\link{iwt1}}, \code{\link{iwt2}}}
\examples{
set.seed(1)
tt <- seq(0, 1, length.out = 60)
y <- matrix(rnorm(30 * 60, sd = 0.5), 30, 60) +
  matrix(exp(-((tt - 0.7) / 0.1)^2), 30, 60, byrow = TRUE)
r <- iwt1(y, n_perm = 200, recycle = FALSE)
plot(r)
g <- autoplot(r, panels = "pvalues", alpha1 = 0.1, title = "One-sample IWT")
r2 <- twt2(y, matrix(rnorm(30 * 60, sd = 0.5), 30, 60), n_perm = 200)
plot(r2, xrange = c(0, 24), xlab = "Hour")
}
