#!/usr/bin/env Rscript
# Package line coverage via covr. Run full test tier before release:
#   INTERACT_SLOW_TESTS=true Rscript scripts/coverage.R

if (!requireNamespace("covr", quietly = TRUE)) {
  stop(
    "Install covr: install.packages(\"covr\")",
    call. = FALSE
  )
}

slow <- identical(Sys.getenv("INTERACT_SLOW_TESTS"), "true")
if (!slow) {
  message(
    "Note: default tests skip slow tier. ",
    "For release coverage run with INTERACT_SLOW_TESTS=true"
  )
}

pkg_root <- if (file.exists("DESCRIPTION")) {
  normalizePath(".")
} else {
  normalizePath(file.path(".."))
}

cov <- covr::package_coverage(
  path = pkg_root,
  type = "tests",
  clean = FALSE
)

pct <- covr::percent_coverage(cov)
message(sprintf("Line coverage: %.1f%%", pct))

by_file <- as.data.frame(cov)
if ("filename" %in% names(by_file)) {
  by_file <- by_file[order(by_file$value), , drop = FALSE]
  names(by_file)[names(by_file) == "value"] <- "coverage"
  print(by_file[, c("filename", "coverage"), drop = FALSE], row.names = FALSE)
}

out_dir <- file.path(pkg_root, "coverage")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
html_path <- file.path(out_dir, "index.html")
covr::report(cov, file = html_path)
message("Wrote ", html_path)

zero <- covr::zero_coverage(cov)
if (length(zero) > 0L) {
  message("\nUncovered lines (first 40 entries):")
  print(head(zero, 40L))
}

invisible(cov)
