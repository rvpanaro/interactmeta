# Information fractions and design tables.

test_that("ifrac() computes IF columns and reference level", {
  res <- ifrac(
    yi = c(0.2, 0.1, 0.3, 0.15),
    vi = c(0.01, 0.04, 0.0225, 0.09),
    trial = c("A", "A", "B", "B"),
    subgroup = c("g1", "g2", "g1", "g2")
  )
  expect_equal(nrow(res), 4L)
  expect_equal(unique(as.character(res$reference)), "g1")
  expect_true("IF_g2" %in% names(res))
  expect_equal(res$IF_g2[[1L]], 0.2)
  expect_true(all(vapply(res[c("yi", "vi", "precision")], is.numeric, logical(1L))))
})

test_that("ifrac() matches escalc data.frame input", {
  skip_if_not_installed("metafor")
  e <- metafor::escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
  res <- ifrac(e, reference = "Lower")
  expect_s3_class(res, "tbl_df")
  expect_true("IF_Higher" %in% names(res))
  expect_equal(unique(as.character(res$reference)), "Lower")
})

test_that("ifrac() on a CAMS fit returns trial-level IF summary", {
  skip_if_not_installed("bayesmeta")
  fit <- fixture_cam_tiny_fit()
  oi <- ifrac(fit)
  expect_s3_class(oi, "interactmeta_overall_ifrac")
  expect_true(length(oi) > 1L)
})