test_that("documented ifrac() example matches named subgroup labels", {
  yi <- c(0.2, 0.1, 0.3, 0.15)
  vi <- c(0.01, 0.04, 0.0225, 0.09)
  trial <- c("A", "A", "B", "B")
  subgroup <- c("g1", "g2", "g1", "g2")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup)

  expect_equal(nrow(res), 4L)
  expect_equal(unique(as.character(res$reference)), "g1")
  expect_true("ifrac_g2" %in% names(res))
  expect_equal(res$ifrac_g2[[1L]], 0.2)
})

test_that("ifrac() returns numeric yi, vi, precision, and ifrac_* columns", {
  res <- ifrac(
    yi = c("0.2", "0.1"),
    vi = c("0.01", "0.04"),
    trial = c("A", "A"),
    subgroup = c("g1", "g2"),
    include_interaction = TRUE
  )
  num_cols <- c("yi", "vi", "precision", grep("^ifrac", names(res), value = TRUE))
  expect_true(all(vapply(res[num_cols], is.numeric, logical(1L))))
  expect_type(res$yi, "double")
})

test_that("ifrac() preserves patient counts in column n", {
  res <- ifrac(edss, reference = "Lower")
  expect_true("n" %in% names(res))
  define <- res[res$study == "DEFINE" & as.character(res$level) == "Lower", , drop = FALSE]
  expect_equal(define$n[[1L]], 336)
  lookup <- interact:::forestplot_trial_n_lookup(res, n_var = "n", contrast_sep = "-")
  expect_equal(lookup[["DEFINE"]], 336 + 73)
})

test_that("ifrac() replaces escalc placeholder rows (yi = 0, vi = 8) with big_var", {
  skip_if_not_installed("metafor")

  dat <- edss
  dat$yi[[1L]] <- 0
  dat$vi[[1L]] <- 8
  e <- metafor::escalc(measure = "IRR", yi = dat$yi, vi = dat$vi, data = dat)

  res <- ifrac(e, reference = "Lower", contrast = NULL)
  row <- res[res$study == dat$study[[1L]] & res$level == dat$level[[1L]], , drop = FALSE]
  expect_equal(row$yi, 0)
  expect_equal(row$vi, 100^2)
})

test_that("ifrac() accepts escalc / data.frame input", {
  skip_if_not_installed("metafor")

  dat <- edss
  e <- metafor::escalc(measure = "IRR", yi = dat$yi, vi = dat$vi, data = dat)

  res_vec <- ifrac(
    yi = dat$yi,
    vi = dat$vi,
    trial = dat$study,
    subgroup = dat$level,
    reference = "Lower"
  )
  res_df <- ifrac(e, reference = "Lower")

  expect_equal(res_df$study, res_vec$study)
  expect_equal(as.character(res_df$level), as.character(res_vec$level))
  expect_equal(res_df$yi, res_vec$yi, tolerance = 1e-12)
  expect_equal(res_df$vi, res_vec$vi, tolerance = 1e-12)
  expect_equal(res_df$ifrac_Higher, res_vec$ifrac_Higher, tolerance = 1e-12)
})

test_that("ifrac() works with different subgroup types and uneven combinations", {
  yi <- c(0.2, 0.1, 0.3)
  vi <- c(0.01, 0.04, 0.09)
  trial <- c("A", "A", "B")
  subgroup_chr <- c("g1", "g2", "g1")

  res_chr <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup_chr, big_var = 1e6, contrast = NULL)
  expect_s3_class(res_chr, "tbl_df")

  # uneven subgroup sets: trial B is missing g2 -> imputed with yi=0 and vi=big_var
  expect_true(any(res_chr$study == "B" & res_chr$level == "g2"))
  row_imp <- res_chr[res_chr$study == "B" & res_chr$level == "g2", , drop = FALSE]
  expect_equal(row_imp$yi, 0)
  expect_equal(row_imp$vi, 1e6)

  # information fractions exist for both subgroup labels (wide output)
  expect_true(all(c("ifrac_g1", "ifrac_g2") %in% names(res_chr)))

  # within each trial, ifrac columns should sum to 1
  sums <- res_chr |>
    dplyr::distinct(.data$study, .data$ifrac_g1, .data$ifrac_g2) |>
    dplyr::mutate(sum_ifrac = .data$ifrac_g1 + .data$ifrac_g2)
  expect_true(all(abs(sums$sum_ifrac - 1) < 1e-12))

  # subgroup as factor should behave identically (column naming may differ in order)
  subgroup_fac <- factor(subgroup_chr, levels = c("g2", "g1"))
  res_fac <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup_fac, big_var = 1e6, contrast = NULL)
  expect_true(all(c("ifrac_g1", "ifrac_g2") %in% names(res_fac)))

  # compare fractions numerically (allow different column order)
  res_chr_frac <- res_chr |>
    dplyr::distinct(.data$study, .data$ifrac_g1, .data$ifrac_g2) |>
    dplyr::arrange(.data$study)
  res_fac_frac <- res_fac |>
    dplyr::distinct(.data$study, .data$ifrac_g1, .data$ifrac_g2) |>
    dplyr::arrange(.data$study)
  expect_equal(res_fac_frac, res_chr_frac)
})

test_that("ifrac() aggregates repeated trial-by-subgroup rows via summed precision", {
  # Two rows for A:g1, one for A:g2
  yi <- c(0.2, 0.25, 0.1)
  vi <- c(0.01, 0.04, 0.04)
  trial <- c("A", "A", "A")
  subgroup <- c("g1", "g1", "g2")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, contrast = NULL)

  # precision sums: g1 prec = 1/0.01 + 1/0.04 = 100 + 25 = 125; g2 prec = 25
  # so fractions should be 125/150 and 25/150
  fr <- res |>
    dplyr::distinct(.data$study, .data$ifrac_g1, .data$ifrac_g2)
  expect_equal(fr$ifrac_g1, 125 / 150, tolerance = 1e-12)
  expect_equal(fr$ifrac_g2, 25 / 150, tolerance = 1e-12)
})

test_that("ifrac() stacks contrast rows with g2-g1 style labels", {
  yi <- c(0.2, 0.1, 0.3)
  vi <- c(0.01, 0.04, 0.09)
  trial <- c("A", "A", "B")
  subgroup <- c("g1", "g2", "g1")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, contrast = c("g1", "g2"), include_interaction = TRUE)

  expect_true(all(c("ifrac_g1") %in% names(res)))
  expect_false(any(grepl("^ifrac2_", names(res))))
  expect_false(any(is.na(res$ifrac_g1)))

  expect_true("g1-g2" %in% res$level)
  expect_equal(nrow(res), 6)

  row_contrast <- res[res$study == "B" & res$level == "g1-g2", , drop = FALSE]
  expect_equal(row_contrast$yi, 0.3 - 0)
  expect_equal(row_contrast$vi, 0.09 + 100^2)

  fr_a <- res[res$study == "A" & res$level == "g1-g2", , drop = FALSE]
  expect_equal(fr_a$precision, 100 + 25, tolerance = 1e-12)
  expect_false("ifrac" %in% names(res))
  expect_equal(fr_a$ifrac_g1, 100 / 125, tolerance = 1e-12)

  row_g1 <- res[res$study == "A" & res$level == "g1", , drop = FALSE]
  expect_equal(row_g1$ifrac_g1, 100 / 125, tolerance = 1e-12)
})

test_that("ifrac() omits overall rows by default", {
  yi <- c(0.2, 0.1, 0.3, 0.15)
  vi <- c(0.01, 0.04, 0.0225, 0.09)
  trial <- c("A", "A", "B", "B")
  subgroup <- c("g1", "g2", "g1", "g2")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, reference = "g1")

  expect_equal(nrow(res), 4L)
  expect_false("Combined" %in% res$level)
})

test_that("ifrac() can append overall rows with include_combined = TRUE", {
  yi <- c(0.2, 0.1, 0.3, 0.15)
  vi <- c(0.01, 0.04, 0.0225, 0.09)
  trial <- c("A", "A", "B", "B")
  subgroup <- c("g1", "g2", "g1", "g2")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, reference = "g1", include_combined = TRUE)

  expect_equal(nrow(res), 6L)
  expect_true(all(c("Combined") %in% res$level))
  expect_true(all(c("ifrac_g2") %in% names(res)))
  expect_false(any(grepl("^ifrac2_", names(res))))

  row_a <- res[res$study == "A" & res$level == "Combined", , drop = FALSE]
  expect_equal(row_a$yi, (0.2 * 100 + 0.1 * 25) / 125, tolerance = 1e-12)
  expect_equal(row_a$vi, 1 / 125, tolerance = 1e-12)
  expect_equal(row_a$precision, 125, tolerance = 1e-12)
  expect_equal(row_a$ifrac_g2, 25 / 125, tolerance = 1e-12)
})

test_that("ifrac() can append overall and interaction rows together", {
  yi <- c(0.2, 0.1, 0.05)
  vi <- c(0.01, 0.04, 0.16)
  trial <- c("A", "A", "A")
  subgroup <- c("g1", "g2", "g3")

  res <- ifrac(
    yi = yi,
    vi = vi,
    trial = trial,
    subgroup = subgroup,
    include_combined = TRUE,
    include_interaction = TRUE
  )

  expect_equal(nrow(res), 6L)
  expect_true(all(c("Combined", "g2-g1", "g3-g1") %in% res$level))
})

test_that("ifrac() omits interaction rows by default", {
  yi <- c(0.2, 0.1, 0.05)
  vi <- c(0.01, 0.04, 0.16)
  trial <- c("A", "A", "A")
  subgroup <- c("g1", "g2", "g3")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup)

  expect_equal(nrow(res), 3)
  expect_false(any(grepl("-", res$level, fixed = TRUE)))
  expect_true(all(c("ifrac_g2", "ifrac_g3") %in% names(res)))
  expect_false(any(grepl("^ifrac2_", names(res))))
  expect_true(all(res$reference == "g1"))
  expect_s3_class(res$level, "factor")
  expect_equal(levels(res$level), c("g1", "g2", "g3"))
  expect_identical(levels(res$level)[[1]], unique(res$reference))
})

test_that("ifrac() puts reference as the first level of subgroup factor", {
  skip_if_not_installed("metafor")

  e <- metafor::escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
  res <- ifrac(e, reference = "Lower")

  expect_s3_class(res$level, "factor")
  expect_equal(levels(res$level)[[1]], "Lower")
  expect_equal(unique(res$reference), "Lower")
})

test_that("edss defaults to Lower reference via factor levels", {
  skip_if_not_installed("metafor")

  expect_equal(levels(edss$subgroup), c("Lower", "Higher"))
  e <- metafor::escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
  res <- ifrac(e)
  expect_equal(unique(as.character(res$reference)), "Lower")
  expect_equal(levels(res$level)[[1]], "Lower")
  expect_true("ifrac_Higher" %in% names(res))
  expect_false("ifrac_Lower" %in% names(res))
})

test_that("ifrac() can append interaction rows with include_interaction = TRUE", {
  yi <- c(0.2, 0.1, 0.05)
  vi <- c(0.01, 0.04, 0.16)
  trial <- c("A", "A", "A")
  subgroup <- c("g1", "g2", "g3")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, include_interaction = TRUE)

  expect_true(all(c("ifrac_g2", "ifrac_g3") %in% names(res)))
  expect_false(any(grepl("^ifrac2_", names(res))))
  expect_false(any(is.na(res$ifrac_g2)))
  expect_false(any(is.na(res$ifrac_g3)))

  expect_true(all(c("g2-g1", "g3-g1") %in% res$level))
  expect_equal(nrow(res), 5)

  row_g2g1 <- res[res$level == "g2-g1", , drop = FALSE]
  expect_equal(row_g2g1$yi, 0.1 - 0.2)
  expect_equal(row_g2g1$vi, 0.04 + 0.01)
  expect_equal(row_g2g1$precision, 25 + 100, tolerance = 1e-12)
  expect_equal(row_g2g1$ifrac_g2, 25 / (100 + 25 + 6.25), tolerance = 1e-12)

  row_g2 <- res[res$level == "g2", , drop = FALSE]
  expect_equal(row_g2$ifrac_g2, 25 / (100 + 25 + 6.25), tolerance = 1e-12)

  row_g3 <- res[res$level == "g3", , drop = FALSE]
  expect_equal(row_g3$ifrac_g3, 6.25 / (100 + 25 + 6.25), tolerance = 1e-12)
})

test_that("ifrac() can use a reference subgroup for contrast=all", {
  yi <- c(0.2, 0.1, 0.05)
  vi <- c(0.01, 0.04, 0.16)
  trial <- c("A", "A", "A")
  subgroup <- c("g1", "g2", "g3")

  res <- ifrac(yi = yi, vi = vi, trial = trial, subgroup = subgroup, reference = 1, include_interaction = TRUE)

  expect_true(all(res$reference == "g1"))
  expect_true(all(c("g2-g1", "g3-g1") %in% res$level))
  expect_false("g2-g3" %in% res$level)
})
