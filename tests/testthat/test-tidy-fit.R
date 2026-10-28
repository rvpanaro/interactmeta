test_that("tidy.interactmeta returns export schema and matches CAMS interaction QC", {
  skip_if_not_installed("metafor")

  fit <- cams(age, reference = "Younger", method = "REML", silent = 2)

  td <- tidy(fit, estimands = c("interaction", "modifiers", "combined"))
  expect_s3_class(td, "tbl_df")
  expect_true(all(c(
    "estimand", "estimate", "std.error", "conf.low", "conf.high",
    "reference", "method", "inference"
  ) %in% names(td)))
  expect_equal(unique(td$method), "CAMS")
  expect_equal(unique(td$reference), "Younger")
  expect_true("interaction" %in% unique(td$estimand))

  ix_row <- td |>
    dplyr::filter(estimand == "interaction", type == "Overall")
  expect_equal(nrow(ix_row), 1L)

  ref <- tidy(
    interactions(fit, type = "Overall", summary = TRUE, exponentiate = FALSE)
  )
  expect_equal(ix_row$estimate[[1L]], ref$estimate[[1L]], tolerance = 1e-8)
  expect_equal(ix_row$conf.low[[1L]], ref$conf.low[[1L]], tolerance = 1e-8)
})

test_that("tidy.interactmeta inherits exponentiate from fit$reporting", {
  skip_if_not_installed("metafor")
  fit <- cams(age, reference = "Younger", method = "REML", silent = 2)
  fit$reporting$exponentiate <- TRUE
  td <- tidy(fit, estimands = "interaction")
  ref <- tidy(
    interactions(fit, type = "Overall", summary = TRUE, exponentiate = TRUE)
  )
  expect_equal(td$estimate[[1L]], ref$estimate[[1L]], tolerance = 1e-6)
})
