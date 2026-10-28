test_that("add_ifrac_design() matches the three-subgroup centered design", {
  dat <- dplyr::tibble(
    trial = "A",
    level = c("g1", "g2", "g3"),
    yi = c(0.2, 0.1, 0.05),
    vi = c(0.01, 0.04, 0.16),
    reference = "g1"
  )

  res <- add_ifrac_design(dat)

  expect_equal(res$precision, c(100, 25, 6.25), tolerance = 1e-10)
  expect_true(is.numeric(res$precision))
  expect_true(all(vapply(res[grep("^ifrac", names(res), value = TRUE)], is.numeric, logical(1L))))
  expect_false("ifrac" %in% names(res))

  expect_equal(res$ifrac_g2, rep(25 / 131.25, 3), tolerance = 1e-10)
  expect_equal(res$ifrac_g3, rep(6.25 / 131.25, 3), tolerance = 1e-10)
  expect_false(any(grepl("^ifrac2_", names(res))))
})

test_that("add_ifrac_design() uses arbitrary subgroup labels", {
  dat <- dplyr::tibble(
    study = rep("T1", 3),
    level = c("low", "medium", "high"),
    yi = c(0, 0, 0),
    vi = c(1, 2, 4),
    reference = "low"
  )

  res <- add_ifrac_design(dat)

  expect_true(all(c("ifrac_medium", "ifrac_high") %in% names(res)))
  expect_false(any(grepl("^ifrac2_", names(res))))
  expect_false(any(grepl("^ifrac_low", names(res))))

  row_med <- res[res$level == "medium", , drop = FALSE]
  expect_equal(row_med$ifrac_medium, (1 / 2) / (1 + 1 / 2 + 1 / 4), tolerance = 1e-12)
})

test_that("add_ifrac2_design() is added only for CAMS fitting", {
  dat <- dplyr::tibble(
    trial = "A",
    level = c("g1", "g2"),
    yi = c(0.2, 0.1),
    vi = c(0.01, 0.04),
    reference = "g1"
  )
  res <- add_ifrac_design(dat)
  fit_dat <- interact:::add_ifrac2_design(res)
  expect_true("ifrac2_g2" %in% names(fit_dat))
  expect_equal(
    fit_dat$ifrac2_g2,
    as.numeric(fit_dat == "g2") - fit_dat$ifrac_g2,
    tolerance = 1e-12
  )
})

test_that("add_ifrac_design() excludes overall rows by default", {
  dat <- dplyr::tibble(
    study = rep("A", 4),
    level = c("g1", "g2", "g3", "Combined"),
    yi = c(0, 0, 0, 0),
    vi = c(0.01, 0.04, 0.16, 0.008),
    reference = "g1"
  )

  res <- add_ifrac_design(dat)
  expect_equal(nrow(res), 3L)
  expect_false("Combined" %in% res$level)
})

test_that("add_ifrac_design() excludes contrast rows by default", {
  dat <- dplyr::tibble(
    study = rep("A", 4),
    level = c("g1", "g2", "g3", "g2-g1"),
    yi = c(0, 0, 0, 0),
    vi = c(0.01, 0.04, 0.16, 0.05),
    reference = "g1"
  )

  res <- add_ifrac_design(dat)
  expect_equal(nrow(res), 3L)
  expect_false("g2-g1" %in% res$level)
})
