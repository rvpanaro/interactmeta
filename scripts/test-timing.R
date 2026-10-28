#!/usr/bin/env Rscript
# Per-file testthat timing baseline. Writes tests/testthat/.timing-baseline.csv

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install testthat.", call. = FALSE)
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}

slow <- identical(Sys.getenv("INTERACT_SLOW_TESTS"), "true")
if (!slow) {
  message("Running fast tier (INTERACT_SLOW_TESTS unset).")
} else {
  message("Running full tier (INTERACT_SLOW_TESTS=true).")
}

pkg_root <- if (file.exists("DESCRIPTION")) "." else ".."
pkgload::load_all(pkg_root, quiet = TRUE)

test_dir <- file.path(pkg_root, "tests", "testthat")
files <- sort(list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE))

rows <- lapply(files, function(f) {
  t0 <- Sys.time()
  err <- NULL
  tryCatch(
    testthat::test_file(f, reporter = "silent"),
    error = function(e) {
      err <<- conditionMessage(e)
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  data.frame(
    file = basename(f),
    seconds = round(elapsed, 2),
    error = err %||% "",
    stringsAsFactors = FALSE
  )
})

tab <- do.call(rbind, rows)
tab <- tab[order(-tab$seconds), , drop = FALSE]
out <- file.path(test_dir, ".timing-baseline.csv")
write.csv(tab, out, row.names = FALSE)
message("Wrote ", out)
print(head(tab, 15L), row.names = FALSE)
message(sprintf("Total: %.1f s", sum(tab$seconds)))
