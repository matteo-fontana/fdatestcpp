# Plot methods: the ggplot is built without warnings for every kind of result,
# the shaded regions match the adjusted p-values, and the abscissa, panels and
# arguments behave as documented.

library(fdatestcpp)

n_checks <- 0L
fails <- character(0)
expect <- function(ok, label) {
  n_checks <<- n_checks + 1L
  if (!isTRUE(ok)) fails <<- c(fails, label)
}
err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
# Build the plot, recording warnings and messages
build <- function(g) {
  msgs <- character(0)
  b <- withCallingHandlers(ggplot2::ggplot_build(g),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  list(b = b, msgs = msgs)
}
geoms <- function(g) vapply(g$layers, function(l) class(l$geom)[1L], "")
# Rectangles of the first panel (they are repeated in every panel)
rect_data <- function(g) {
  i <- which(geoms(g) == "GeomRect")
  if (!length(i)) return(NULL)
  d <- ggplot2::layer_data(g, i)
  d[d$PANEL == 1L, ]
}

# Shaded regions on a hand-made result
fake <- function(adj, cls = "fos", ...) {
  p <- length(adj)
  obj <- list(data = matrix(seq_len(3 * p), 3, p), mu = 0, unadjusted_pvalues = adj / 2,
              adjusted_pvalues = adj, correction_method = "IWT", ...)
  if (cls == "fts") obj$group_labels <- c(1, 1, 2)
  structure(obj, class = c(paste0("fdatestcpp_", cls), cls))
}
g <- autoplot(fake(c(0.5, 0.03, 0.005, 0.004, 0.2)))
rd <- rect_data(g)
expect(isTRUE(all.equal(rd$xmin, c(0.125, 0.375))) && isTRUE(all.equal(rd$xmax, c(0.375, 0.875))),
       "rectangles span the cells of the significant points")
expect(identical(rd$fill, c("grey89", "grey76")), "two shades, darker below alpha2")
rd <- rect_data(autoplot(fake(c(0.001, 0.5, 0.5, 0.5, 0.001))))
expect(isTRUE(all.equal(rd$xmin, c(0, 0.875))) && isTRUE(all.equal(rd$xmax, c(0.125, 1))),
       "rectangles cut at the ends of the domain")
expect(is.null(rect_data(autoplot(fake(rep(0.5, 5))))), "nothing significant: no rectangles")
rd <- rect_data(autoplot(fake(c(0.04, 0.04, 0.001)), alpha1 = 0.01, alpha2 = 0.05))
expect(nrow(rd) == 2L && isTRUE(all.equal(rd$xmax, c(0.75, 1))), "alpha1 < alpha2 swapped")
rd <- rect_data(autoplot(fake(c(0.04, 0.001, 0.5)), alpha1 = 0.05, alpha2 = 0.05))
expect(nrow(rd) == 1L && isTRUE(all.equal(c(rd$xmin, rd$xmax), c(0, 0.75))), "equal levels: one shade")
g <- autoplot(fake(c(0.04999999995, 0.001, 0.5)), alpha1 = 0.05, alpha2 = 0.0499999999)
expect(nrow(rect_data(g)) == 2L, "levels that print alike")
rd <- rect_data(autoplot(fake(c(NA, 0.001, NaN, 0.001))))
expect(nrow(rd) == 2L, "NA p-values are not shaded")
rd <- rect_data(autoplot(fake(c(0.5, 0.001, 0.5)), xrange = c(10, 20)))
expect(isTRUE(all.equal(c(rd$xmin, rd$xmax), c(12.5, 17.5))), "xrange")
rd <- rect_data(autoplot(fake(c(0.5, 0.001, 0.5), abscissa = c(0, 1, 4))))
expect(isTRUE(all.equal(c(rd$xmin, rd$xmax), c(0.5, 2.5))), "uneven abscissa")

# Real results: every wrapper, matrix and fd data, all panel choices
set.seed(3)
tt <- seq(0, 1, length.out = 40)
y <- matrix(rnorm(20 * 40, sd = 0.4), 20, 40) +
  matrix(exp(-((tt - 0.7) / 0.1)^2), 20, 40, byrow = TRUE)
z <- matrix(rnorm(20 * 40, sd = 0.4), 20, 40)
basis <- fda::create.bspline.basis(c(2, 6), nbasis = 12)
y_fd <- fda::smooth.basis(2 + 4 * tt, t(y), basis)$fd
z_fd <- fda::smooth.basis(2 + 4 * tt, t(z), basis)$fd
B <- 100L
results <- list(
  iwt1 = iwt1(y, n_perm = B, recycle = FALSE), twt1 = twt1(y, n_perm = B),
  global1 = global1(y, n_perm = B), iwt2 = iwt2(y, z, n_perm = B, recycle = FALSE),
  twt2 = twt2(y, z, n_perm = B), global2 = global2(y, z, n_perm = B),
  paired = iwt2(y, z, n_perm = B, paired = TRUE), fd1 = iwt1(y_fd, n_perm = B),
  fd2 = twt2(y_fd, z_fd, n_perm = B),
  fd_mu = iwt1(y_fd, mu = fda::fd(matrix(0.2, 12, 1), basis), n_perm = B),
  mu2 = iwt2(y, z, mu = 0.3, n_perm = B),
  nan = iwt1(cbind(0, y[, 1:10], 0), standardize = TRUE, n_perm = B),
  p1 = iwt2(y[, 28, drop = FALSE], z[, 28, drop = FALSE], n_perm = B),
  p1_one = twt1(y[, 28, drop = FALSE], n_perm = B))
one_sample <- c("iwt1", "twt1", "global1", "fd1", "fd_mu", "nan", "p1_one")
panel_sets <- list(c("data", "pvalues"), "data", "pvalues", c("pvalues", "data"))
for (nm in names(results)) {
  r <- results[[nm]]
  expect(identical(class(r), if (nm %in% one_sample) c("fdatestcpp_fos", "fos") else
    c("fdatestcpp_fts", "fts")), paste(nm, "class"))
  for (ps in panel_sets) for (un in c(TRUE, FALSE)) {
    g <- autoplot(r, panels = ps, unadjusted = un, means = un)
    lab <- sprintf("%s panels=%s unadjusted=%s", nm, paste(ps, collapse = "+"), un)
    expect(inherits(g, "ggplot"), paste(lab, "is ggplot"))
    res <- build(g)
    expect(!length(res$msgs), paste(lab, "no warnings:", paste(res$msgs, collapse = "; ")))
    expect(nrow(res$b$layout$layout) == length(ps), paste(lab, "number of panels"))
  }
}
# fd data: drawn on the evaluation grid
g <- autoplot(results$fd1, panels = "data", means = FALSE)
expect(isTRUE(all.equal(sort(unique(ggplot2::layer_data(g, length(g$layers))$x)),
                        results$fd1$abscissa)), "fd abscissa")
expect(isTRUE(all.equal(range(results$fd1$abscissa), c(2, 6))), "fd range")
# Global results: the whole domain is shaded when the global test rejects
r <- results$global1
rd <- rect_data(autoplot(r))
if (r$global_pvalue < 0.05) {
  expect(nrow(rd) == 1L && isTRUE(all.equal(c(rd$xmin, rd$xmax), c(0, 1))), "global shading")
}

# plot() draws and returns the ggplot invisibly
grDevices::pdf(NULL)
vis <- withVisible(plot(results$iwt2, title = "test", alpha1 = 0.1))
grDevices::dev.off()
expect(!vis$visible && inherits(vis$value, "ggplot"), "plot() returns the ggplot invisibly")
expect(identical(vis$value$labels$title, "test") ||
         identical(ggplot2::get_labs(vis$value)$title, "test"), "title passed through")

# Coarse grids: p-values drawn as steps over the shaded cells; fine grids: lines
pv_layer <- function(g) {
  i <- which(vapply(g$layers, function(l) !is.null(l$mapping$linetype), NA))
  list(geom = class(g$layers[[i]]$geom)[1L], data = ggplot2::layer_data(g, i))
}
L <- pv_layer(autoplot(fake(c(0.5, 0.03, 0.005, 0.004, 0.2)), panels = "pvalues"))
adj <- L$data[L$data$linetype == "solid", ]
expect(L$geom == "GeomPath", "coarse grid: steps")
expect(isTRUE(all.equal(adj$x, c(0, 0.125, 0.125, 0.375, 0.375, 0.625, 0.625, 0.875, 0.875, 1))) &&
         isTRUE(all.equal(adj$y, rep(c(0.5, 0.03, 0.005, 0.004, 0.2), each = 2))),
       "steps change level at the cell edges")
expect(pv_layer(autoplot(fake(c(0.5, 0.03, 0.005)), steps = FALSE))$geom == "GeomLine",
       "steps = FALSE")
expect(pv_layer(autoplot(fake(rep(0.5, 51))))$geom == "GeomLine", "51 points: lines")
expect(pv_layer(autoplot(fake(rep(0.5, 51)), steps = TRUE))$geom == "GeomPath", "steps = TRUE")
res <- build(autoplot(fake(c(0.5, NA, 0.001, 0.2))))
expect(!length(res$msgs), "steps with a missing p-value: no warnings")
expect(err(autoplot(fake(c(0.5, 0.2)), steps = NA)), "steps = NA")

# Two-sample colours, interleaved drawing order and the H0 reference for mu != 0
g <- autoplot(results$iwt2, panels = "data")
ld <- lapply(seq_along(g$layers), function(i) ggplot2::layer_data(g, i))
means_layer <- ld[[length(ld)]]
expect(setequal(unique(means_layer$colour), c("#E69F00", "#0072B2")), "Okabe-Ito sample colours")
curves_layer <- ld[[which(geoms(g) == "GeomLine")[1L]]]
first_groups <- curves_layer$colour[match(sort(unique(curves_layer$group))[1:4], curves_layer$group)]
expect(length(unique(first_groups)) == 2L, "samples interleaved")
g <- autoplot(results$mu2, panels = "data")
dashed <- vapply(g$layers, function(l) identical(l$aes_params$linetype, "dashed"), NA)
expect(sum(dashed) == 1L, "two-sample mu != 0: dashed reference line")
ref <- ggplot2::layer_data(g, which(dashed))
r <- results$mu2
m2 <- colMeans(r$data[r$group_labels == 2, ])
expect(isTRUE(all.equal(ref$y, m2 + 0.3)), "reference = mean of sample 2 + mu")
g <- autoplot(results$iwt2, panels = "data")
expect(!any(vapply(g$layers, function(l) identical(l$aes_params$linetype, "dashed"), NA)),
       "two-sample mu = 0: no reference line")
caption <- function(g) {
  cap <- g$labels$caption
  if (is.null(cap) && exists("get_labs", asNamespace("ggplot2"))) cap <- ggplot2::get_labs(g)$caption
  cap
}
expect(grepl("mean of sample 2 + mu", caption(autoplot(results$mu2)), fixed = TRUE),
       "caption explains the dashed line")
expect(is.null(caption(autoplot(results$iwt2))), "no caption when mu = 0")
# a two-colour scale recolours samples, means and reference line
g <- suppressMessages(autoplot(results$mu2, panels = "data") +
  ggplot2::scale_colour_manual(values = c("1" = "red", "2" = "darkgreen")))
res <- build(g)
dashed <- which(vapply(g$layers, function(l) identical(l$aes_params$linetype, "dashed"), NA))
expect(!length(res$msgs) && identical(unique(ggplot2::layer_data(g, dashed)$colour), "red"),
       "colour override recolours the reference line")
# sample 1 (lighter orange) curves more opaque than sample 2
g <- autoplot(results$iwt2, panels = "data")
ld <- ggplot2::layer_data(g, which(geoms(g) == "GeomLine")[1L])
expect(length(unique(ld$alpha[ld$colour == "#E69F00"])) == 1L &&
         unique(ld$alpha[ld$colour == "#E69F00"]) > unique(ld$alpha[ld$colour == "#0072B2"]),
       "curve opacity balanced between the samples")
lab <- fdatestcpp:::level_labels(0.05, 0.0499999999)
expect(length(unique(lab)) == 2L, "levels that print alike get distinct labels")
expect(identical(unname(fdatestcpp:::level_labels(0.05, 0.01)), c("< 0.01", "< 0.05")),
       "usual labels")

# xrange moving fd data off their grid warns; matching or matrix data do not
warns <- function(expr) {
  w <- NULL
  withCallingHandlers(expr, warning = function(c) { w <<- conditionMessage(c); invokeRestart("muffleWarning") })
  w
}
w <- warns(autoplot(results$fd1, xrange = c(0, 1)))
expect(is.character(w) && grepl("xrange", w), "xrange on fd data warns")
expect(is.null(warns(autoplot(results$fd1, xrange = range(results$fd1$abscissa)))),
       "xrange equal to the fd grid: no warning")
expect(is.null(warns(autoplot(results$iwt1, xrange = c(0, 24)))), "xrange on matrix data: no warning")

# Extra arguments are reported, not silently dropped
w <- tryCatch(autoplot(results$iwt1, main = "x"), warning = function(w) conditionMessage(w))
expect(is.character(w) && grepl("main", w), "warning for ignored arguments")
grDevices::pdf(NULL)
w <- tryCatch(plot(results$iwt2, lwd = 2), warning = function(w) conditionMessage(w))
grDevices::dev.off()
expect(is.character(w) && grepl("lwd", w), "plot() passes extra arguments to the check")

# Panels and labels
expect(nrow(build(autoplot(results$iwt1, panels = c("data", "data")))$b$layout$layout) == 1L,
       "repeated panel")
expect(!err(autoplot(results$iwt1, ylabel = "p-value")), "ylabel equal to the p-value strip")
expect(err(autoplot(results$iwt1, ylabel = NULL)), "NULL ylabel")
expect(err(autoplot(results$iwt1, ylabel = c("a", "b"))), "two ylabels")

# Argument checks
r <- results$iwt1
expect(err(autoplot(r, alpha1 = 0)), "alpha1 = 0")
expect(err(autoplot(r, alpha2 = c(0.01, 0.05))), "alpha2 length")
expect(err(autoplot(r, alpha1 = 1.5)), "alpha1 > 1")
expect(err(autoplot(r, xrange = c(1, 0))), "decreasing xrange")
expect(err(autoplot(r, xrange = c(0, NA))), "NA xrange")
expect(err(autoplot(r, panels = "matrix")), "unknown panel")
expect(err(autoplot(r, unadjusted = NA)), "NA flag")
expect(err(autoplot(r, linewidth = -1)), "negative linewidth")
bad <- r
bad$data <- bad$data[, -1]
expect(err(autoplot(bad)), "data and p-values of different lengths")

cat(sprintf("test-plot: %d checks, %d failed\n", n_checks, length(fails)))
if (length(fails)) stop("Failed:\n", paste(fails, collapse = "\n"))
