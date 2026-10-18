# ggplot2 graphics for the results of iwt1(), twt1(), global1(), iwt2(),
# twt2() and global2(), after the IWTgg plotting scheme (M. Fontana, 2018):
# the curves with the part of the domain selected by the test shaded in grey
# (adjusted p-value below alpha1, darker below alpha2), and the adjusted
# (solid) and unadjusted (dashed) p-value functions with both levels drawn as
# horizontal lines. The panels are facets of one ggplot over a common abscissa,
# so the result can be modified, combined and saved as any ggplot.

autoplot.fdatestcpp_fos <- function(object, panels = c("data", "pvalues"), alpha1 = 0.05,
                                    alpha2 = 0.01, xrange = NULL, unadjusted = TRUE,
                                    steps = NULL, means = TRUE, title = NULL,
                                    xlab = "Domain", ylabel = "Functional data",
                                    linewidth = 0.5, ...) {
  warn_ignored(match.call(expand.dots = FALSE)$...)
  gg_result(object, panels, alpha1, alpha2, xrange, unadjusted, steps, means, title, xlab,
            ylabel, linewidth)
}

autoplot.fdatestcpp_fts <- function(object, panels = c("data", "pvalues"), alpha1 = 0.05,
                                    alpha2 = 0.01, xrange = NULL, unadjusted = TRUE,
                                    steps = NULL, means = TRUE, title = NULL,
                                    xlab = "Domain", ylabel = "Functional data",
                                    linewidth = 0.5, ...) {
  warn_ignored(match.call(expand.dots = FALSE)$...)
  gg_result(object, panels, alpha1, alpha2, xrange, unadjusted, steps, means, title, xlab,
            ylabel, linewidth)
}

plot.fdatestcpp_fos <- function(x, ...) {
  g <- autoplot.fdatestcpp_fos(x, ...)
  print(g)
  invisible(g)
}

plot.fdatestcpp_fts <- function(x, ...) {
  g <- autoplot.fdatestcpp_fts(x, ...)
  print(g)
  invisible(g)
}

# Colours of the two samples (Okabe-Ito orange and blue: distinct for
# colour-blind readers and in greyscale) and of mu in one-sample plots.
sample_colours <- c("1" = "#E69F00", "2" = "#0072B2")
mu_colour <- "#0072B2"

warn_ignored <- function(dots) {
  if (!length(dots)) return(invisible())
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  nm[nm == ""] <- "<unnamed>"
  warning("Arguments ignored by autoplot(): ", paste(nm, collapse = ", "), ".", call. = FALSE)
}

check_level <- function(a, name) {
  if (!is.numeric(a) || length(a) != 1L || is.na(a) || a <= 0 || a > 1) {
    stop(sprintf("`%s` must be a number in (0, 1].", name), call. = FALSE)
  }
  as.double(a)
}

check_label <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single string.", name), call. = FALSE)
  }
  x
}

# Grid on which the p-values are given: `xrange` as in fdatest (p equally
# spaced points), else the abscissa of fd data, else fdatest's [0, 1]. An
# `xrange` that moves the points of fd data off their evaluation grid warns.
plot_abscissa <- function(object, xrange, p) {
  grid <- object$abscissa
  if (!is.null(xrange)) {
    if (!is.numeric(xrange) || length(xrange) != 2L || !all(is.finite(xrange)) ||
        xrange[1L] > xrange[2L] || (p > 1L && xrange[1L] == xrange[2L])) {
      stop("`xrange` must be two increasing finite numbers.", call. = FALSE)
    }
    x <- seq(xrange[1L], xrange[2L], length.out = p)
    if (!is.null(grid) && !isTRUE(all.equal(x, as.double(grid)))) {
      warning(paste0("`xrange` replaces the evaluation grid of the fd data: the p-values are ",
                     "drawn at ", p, " equally spaced points of [", format(xrange[1L]), ", ",
                     format(xrange[2L]), "], not at `abscissa`. To relabel the axis, transform ",
                     "`abscissa` instead (e.g. r$abscissa <- 2000 + 10 * r$abscissa)."),
              call. = FALSE)
    }
    return(x)
  }
  if (!is.null(grid)) return(as.double(grid))
  seq(0, 1, length.out = p)
}

# Cells of the grid points: from the midpoint to the previous point to the
# midpoint to the next one, cut at the first and last points.
cells <- function(x) {
  p <- length(x)
  if (p == 1L) return(list(lo = x - 0.5, hi = x + 0.5))
  mid <- (x[-1L] + x[-p]) / 2
  list(lo = c(x[1L], mid), hi = c(mid, x[p]))
}

# A p-value function as steps that are constant over each cell, the same
# cells as the shading, so the drawn value at any abscissa is the p-value of
# the grid point whose cell contains it. A missing p-value breaks the path.
stairs <- function(x, y) {
  cl <- cells(x)
  data.frame(x = as.vector(rbind(cl$lo, cl$hi)), y = rep(y, each = 2L))
}

# Labels of the shading levels, keyed by level: "2" below alpha2 (darker),
# "1" below alpha1. Equal levels collapse into "2"; levels that differ get
# as many digits as it takes to print them differently.
level_labels <- function(alpha1, alpha2) {
  if (alpha1 == alpha2) return(c("2" = paste("<", format(alpha2, digits = 3))))
  for (d in 3:17) {
    a2 <- format(alpha2, digits = d)
    a1 <- format(alpha1, digits = d)
    if (a1 != a2) break
  }
  c("2" = paste("<", a2), "1" = paste("<", a1))
}

# Shaded regions: runs of consecutive grid points at the same significance
# level, each point covering its cell.
sig_regions <- function(x, pv, alpha1, alpha2) {
  lvl <- ifelse(!is.na(pv) & pv < alpha2, 2L, ifelse(!is.na(pv) & pv < alpha1, 1L, 0L))
  cl <- cells(x)
  lo <- cl$lo
  hi <- cl$hi
  runs <- rle(lvl)
  end <- cumsum(runs$lengths)
  start <- end - runs$lengths + 1L
  keep <- runs$values > 0L
  data.frame(xmin = lo[start[keep]], xmax = hi[end[keep]],
             level = factor(as.character(runs$values[keep]), levels = c("2", "1")))
}

# Curves as lines; with a single grid point, as points (a mapped linetype
# becomes the shape).
trace_layer <- function(single, data, mapping, linewidth, dashed = FALSE, ...) {
  if (single) {
    if (!is.null(mapping$linetype)) {
      mapping$shape <- mapping$linetype
      mapping$linetype <- NULL
    }
    ggplot2::geom_point(data = data, mapping = mapping, size = 1 + 2 * linewidth, ...)
  } else if (dashed) {
    ggplot2::geom_line(data = data, mapping = mapping, linewidth = linewidth,
                       linetype = "dashed", ...)
  } else {
    ggplot2::geom_line(data = data, mapping = mapping, linewidth = linewidth, ...)
  }
}

gg_result <- function(object, panels, alpha1, alpha2, xrange, unadjusted, steps, means,
                      title, xlab, ylabel, linewidth) {
  panels <- unique(match.arg(panels, c("data", "pvalues"), several.ok = TRUE))
  alpha1 <- check_level(alpha1, "alpha1")
  alpha2 <- check_level(alpha2, "alpha2")
  if (alpha1 < alpha2) {
    tmp <- alpha1
    alpha1 <- alpha2
    alpha2 <- tmp
  }
  check_flag(unadjusted)
  if (!is.null(steps)) check_flag(steps)
  check_flag(means)
  check_label(ylabel, "ylabel")
  if (!is.numeric(linewidth) || length(linewidth) != 1L || !is.finite(linewidth) ||
      linewidth <= 0) {
    stop("`linewidth` must be a positive number.", call. = FALSE)
  }
  data <- object$data
  adj <- as.double(object$adjusted_pvalues)
  p <- length(adj)
  if (!is.matrix(data) || ncol(data) != p || p < 1L) {
    stop("`object` must be a result of iwt1(), twt1(), global1(), iwt2(), twt2() or global2().",
         call. = FALSE)
  }
  x <- plot_abscissa(object, xrange, p)
  single <- p == 1L
  # coarse grids: p-values as steps over the cells, so that the drawn curve
  # crosses a level exactly where the shading starts
  if (is.null(steps)) steps <- p <= 50L
  method <- object$correction_method
  global <- identical(method, "Global")
  shade_title <- if (global) "Global p-value" else paste(method, "adjusted p-value")
  pv_labels <- c(adjusted = if (global) "global" else paste(method, "adjusted"),
                 unadjusted = "unadjusted")
  panel_of <- function(panel) factor(panel, levels = panels)
  guide_list <- list()

  g <- ggplot2::ggplot()
  rects <- sig_regions(x, adj, alpha1, alpha2)
  if (nrow(rects)) {
    brk <- intersect(c("2", "1"), as.character(rects$level))
    g <- g +
      ggplot2::geom_rect(data = rects,
                         ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf,
                                      ymax = Inf, fill = .data$level), colour = NA) +
      ggplot2::scale_fill_manual(values = c("2" = "grey76", "1" = "grey89"), breaks = brk,
                                 labels = unname(level_labels(alpha1, alpha2)[brk]),
                                 name = shade_title)
    guide_list$fill <- ggplot2::guide_legend(order = 1)
  }

  if ("data" %in% panels) {
    n <- nrow(data)
    curve_alpha <- max(0.05, min(0.35, 15 / n))
    halo <- 3 * linewidth + 1
    if (inherits(object, "fts")) {
      grp <- factor(object$group_labels, levels = c(1, 2), labels = c("1", "2"))
      # interleave the two samples, so that neither is drawn over the other
      key <- stats::ave(seq_len(n), grp, FUN = function(i) seq_along(i) / length(i))
      id <- integer(n)
      id[order(key)] <- seq_len(n)
      curves <- data.frame(x = rep(x, times = n), y = as.vector(t(data)),
                           id = rep(id, each = p), group = rep(grp, each = p),
                           panel = panel_of("data"))
      # the orange of sample 1 is lighter than the blue of sample 2: more
      # opacity for its curves gives both clouds the same contrast on white
      g <- g + trace_layer(single, curves,
                           ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                        colour = .data$group, alpha = .data$group),
                           linewidth, show.legend = if (means) FALSE else NA) +
        ggplot2::scale_alpha_manual(values = c("1" = min(1, 1.75 * curve_alpha), "2" = curve_alpha),
                                    guide = "none")
      if (means) {
        m1 <- colMeans(data[grp == "1", , drop = FALSE])
        m2 <- colMeans(data[grp == "2", , drop = FALSE])
        avg <- data.frame(x = rep(x, 2L), y = c(m1, m2), group = factor(rep(c("1", "2"), each = p)),
                          panel = panel_of("data"))
        mu <- rep_len(as.double(object$mu), p)
        if (any(mu != 0)) {
          # under H0 the mean of sample 1 is the mean of sample 2 plus mu: drawn
          # dashed in sample 1's colour (so a user colour scale recolours it
          # too) and explained in the caption rather than in the legend
          ref <- data.frame(x = x, y = m2 + mu, group = factor("1", levels = c("1", "2")),
                            panel = panel_of("data"))
          g <- g +
            trace_layer(single, ref, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group),
                        2 * linewidth, dashed = TRUE, show.legend = FALSE) +
            ggplot2::labs(caption = "Dashed: mean of sample 2 + mu, the mean of sample 1 under H0")
        }
        g <- g +
          trace_layer(single, avg, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
                      halo, colour = "white") +
          trace_layer(single, avg, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group),
                      3 * linewidth)
      }
      g <- g + ggplot2::scale_colour_manual(values = sample_colours, name = "Sample")
      guide_list$colour <- ggplot2::guide_legend(
        order = 2, override.aes = if (single) list(alpha = 1) else
          list(alpha = 1, linewidth = 1.5 * linewidth + 0.5))
    } else {
      curves <- data.frame(x = rep(x, times = n), y = as.vector(t(data)),
                           id = rep(seq_len(n), each = p), panel = panel_of("data"))
      g <- g + trace_layer(single, curves,
                           ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                           linewidth, colour = "grey20", alpha = curve_alpha)
      if (means) {
        avg <- data.frame(x = x, y = colMeans(data), line = "mean", panel = panel_of("data"))
        ref <- data.frame(x = x, y = rep_len(as.double(object$mu), p), line = "mu",
                          panel = panel_of("data"))
        g <- g +
          trace_layer(single, avg, ggplot2::aes(x = .data$x, y = .data$y), halo,
                      colour = "white") +
          trace_layer(single, avg, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$line),
                      3 * linewidth) +
          trace_layer(single, ref, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$line),
                      2 * linewidth, dashed = TRUE) +
          ggplot2::scale_colour_manual(values = c(mean = "black", mu = mu_colour),
                                       breaks = c("mean", "mu"),
                                       labels = c("sample mean", "mu (H0)"), name = NULL)
        guide_list$colour <- ggplot2::guide_legend(order = 2)
      }
    }
  }

  if ("pvalues" %in% panels) {
    type <- if (unadjusted) c("adjusted", "unadjusted") else "adjusted"
    vals <- list(adjusted = adj, unadjusted = as.double(object$unadjusted_pvalues))[type]
    pv <- do.call(rbind, lapply(type, function(tp) {
      d <- if (steps && !single) stairs(x, vals[[tp]]) else data.frame(x = x, y = vals[[tp]])
      d$type <- factor(tp, levels = c("adjusted", "unadjusted"))
      d
    }))
    pv$panel <- panel_of("pvalues")
    hl <- data.frame(y = unique(c(alpha1, alpha2)), panel = panel_of("pvalues"))
    g <- g +
      ggplot2::geom_hline(data = hl, ggplot2::aes(yintercept = .data$y), colour = "grey20",
                          linewidth = linewidth, linetype = "dotted") +
      ggplot2::geom_blank(data = data.frame(x = x[1L], y = c(0, 1),
                                            panel = panel_of("pvalues")),
                          ggplot2::aes(x = .data$x, y = .data$y)) +
      (if (steps && !single) {
        ggplot2::geom_path(data = pv, ggplot2::aes(x = .data$x, y = .data$y,
                                                   linetype = .data$type),
                           linewidth = 2 * linewidth, na.rm = TRUE)
      } else {
        trace_layer(single, pv, ggplot2::aes(x = .data$x, y = .data$y, linetype = .data$type),
                    2 * linewidth, na.rm = TRUE)
      })
    if (single) {
      g <- g + ggplot2::scale_shape_manual(values = c(adjusted = 16, unadjusted = 1),
                                           breaks = type, labels = unname(pv_labels[type]),
                                           name = NULL)
      guide_list$shape <- if (unadjusted) ggplot2::guide_legend(order = 3) else "none"
    } else {
      g <- g + ggplot2::scale_linetype_manual(values = c(adjusted = "solid", unadjusted = "dashed"),
                                              breaks = type, labels = unname(pv_labels[type]),
                                              name = NULL)
      guide_list$linetype <- if (unadjusted) ggplot2::guide_legend(order = 3) else "none"
    }
  }
  if (single) g <- g + ggplot2::scale_x_continuous(breaks = x)

  g +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", switch = "y",
                        labeller = ggplot2::as_labeller(c(data = ylabel, pvalues = "p-value"))) +
    do.call(ggplot2::guides, guide_list) +
    ggplot2::labs(title = title, x = xlab, y = NULL) +
    theme_fdatestcpp()
}

# IWTgg's theme: minimal, no strip background, grey grid; the facet strips
# on the left serve as y-axis titles, and the legends are stacked below the
# plot so that they fit the width of a single figure.
theme_fdatestcpp <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 90, size = ggplot2::rel(1.25)),
                   # grid drawn over the layers in a translucent dark grey, so
                   # that it stays visible on the shaded cells
                   panel.ontop = TRUE,
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "#00000024"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.location = "plot",
                   legend.box = "vertical",
                   legend.box.just = "left",
                   legend.spacing.y = ggplot2::unit(2, "pt"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.key.width = ggplot2::unit(2.2, "lines"))
}
