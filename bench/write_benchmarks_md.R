# Writes benchmarks.md from bench/results.csv (and bench/stages.csv if present).
# Run from the fdatestcpp folder.

suppressPackageStartupMessages(library(data.table))

res <- fread(file.path("bench", "results.csv"))
res[, fun := factor(fun, c("iwt1", "iwt2_paired", "twt2_paired"))]
setorder(res, fun, agg, p, B)
lab <- c(iwt1 = "iwt1", iwt2_paired = "iwt2 (paired)", twt2_paired = "twt2 (paired)")
fmt_s <- function(x) ifelse(x >= 10, sprintf("%.0f", x), ifelse(x >= 1, sprintf("%.1f", x),
                     ifelse(x >= 0.01, sprintf("%.3f", x), sprintf("%.4f", x))))
tab <- res[, .(Test = lab[as.character(fun)], Aggregation = agg, p, B,
               `fdatest (s)` = fmt_s(fdatest_s), `fdatestcpp (s)` = fmt_s(fdatestcpp_s),
               `Speed-up` = paste0(formatC(round(speedup), format = "d", big.mark = ","), "x"),
               `Identical p-values` = ifelse(identical, "yes", "NO"))]

md_table <- function(d, left = c("Test", "Aggregation", "Stage", "Identical p-values")) {
  head <- paste0("| ", paste(names(d), collapse = " | "), " |")
  align <- paste0("|", paste(ifelse(names(d) %in% left, ":---", "---:"), collapse = "|"), "|")
  rows <- apply(d, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |"))
  c(head, align, rows)
}

cpu <- Sys.getenv("PROCESSOR_IDENTIFIER")
fd <- packageDescription("fdatest")
lines <- c(
  "# Benchmarks: fdatest vs fdatestcpp",
  "",
  sprintf(paste("Run on %s: %s, Windows 11, Intel Core i7-1360P laptop (%s), single thread",
                "pinned to the performance cores (`start /affinity FF`), no mirai daemons."),
          format(Sys.Date()), R.version.string, if (nzchar(cpu)) cpu else R.version$platform),
  sprintf("fdatest %s from %s, commit %s; fdatestcpp %s.", fd$Version,
          fd$Repository, substr(fd$RemoteSha, 1, 7), packageVersion("fdatestcpp")),
  "",
  paste("End to end, same data and seed: n = 240 curves (pairs for iwt2/twt2), p grid",
        "points, B resamples, `recycle = FALSE`; simulated curves (FAR(1) noise, a bump",
        "in the mean on part of the domain). fdatest is timed once per cell, or as the median of",
        "three runs when a run takes under 5 s; fdatestcpp is the median over repeated",
        "calls. \"Identical\" means `identical()` unadjusted and adjusted p-values and,",
        "for the IWT, the whole p-value matrix."),
  "",
  md_table(tab),
  "",
  sprintf("Speed-up range: %s to %s; median %s.",
          paste0(formatC(round(min(res$speedup)), big.mark = ",", format = "d"), "x"),
          paste0(formatC(round(max(res$speedup)), big.mark = ",", format = "d"), "x"),
          paste0(formatC(round(median(res$speedup)), big.mark = ",", format = "d"), "x")),
  sprintf("Total time for the grid: fdatest %.0f s, fdatestcpp %.1f s.",
          sum(res$fdatest_s), sum(res$fdatestcpp_s)),
  "",
  paste("The IWT gains grow with p because fdatest's adjustment loops over all",
        "p(p+1)/2 intervals in R (`rowSums()` on a copied submatrix for the integral,",
        "`apply(..., 1, max)` for the max), O(B p^3) work with a large constant. The",
        "TWT gains are smaller because fdatest's TWT adjustment is already vectorized",
        "over resamples; there the port's time is mostly the paired resampling (see",
        "below), which fdatest does with an R loop over pairs.")
)

fd_csv <- file.path("bench", "fd_results.csv")
if (file.exists(fd_csv)) {
  fdr <- fread(fd_csv)
  fdr[, fun := factor(fun, c("iwt1", "iwt2_paired", "twt2_paired"))]
  setorder(fdr, fun, agg, B)
  fdtab <- fdr[, .(Test = lab[as.character(fun)], Aggregation = agg, B,
                   `fdatest (s)` = fmt_s(fdatest_s), `fdatestcpp (s)` = fmt_s(fdatestcpp_s),
                   `same curves as matrix (s)` = fmt_s(matrix_input_s),
                   `Speed-up` = paste0(formatC(round(speedup), format = "d", big.mark = ","), "x"),
                   `Identical p-values` = ifelse(identical, "yes", "NO"))]
  lines <- c(lines, "",
    "## fd input",
    "",
    paste0("n = ", fdr$n[1], " fd curves (B-spline basis with 20 functions) evaluated on ",
           "fdatest's default grid (", fdr$p[1], " points), `recycle = FALSE`, same seed. ",
           "The \"same curves as matrix\" column times fdatestcpp on the evaluated curves, so ",
           "the difference to the fd column is the cost of the evaluation. \"Identical\" ",
           "also covers the evaluated `data` field."),
    "",
    md_table(fdtab))
}

stages_csv <- file.path("bench", "stages.csv")
if (file.exists(stages_csv)) {
  st <- fread(stages_csv)
  st[, seconds := fmt_s(seconds)]
  lines <- c(lines, "",
    "## Where fdatestcpp spends its time",
    "",
    paste("Single calls, n = 240, p = 200, B = 5000 (median of repeated calls).",
          "Resampling is the front end (`perm_stats1()` / `perm_stats2()`); the",
          "adjustments run on its output."),
    "",
    md_table(st[, .(Stage = stage, Seconds = seconds)]))
}
writeLines(lines, "benchmarks.md")
cat("wrote benchmarks.md\n")
