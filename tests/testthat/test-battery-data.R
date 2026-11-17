# Packaged data and escalc / HR helpers.

test_that("exported analysis datasets are present", {
  items <- utils::data(package = "interactmeta")$results[, "Item"]
  expect_true(all(c(
    "edss", "age", "il6", "steroids", "respire_age", "respire_sex"
  ) %in% items))
})

test_that("edss and age support ifrac() without error", {
  expect_gt(nrow(ifrac(edss, reference = "Lower")), 4L)
  expect_gt(nrow(ifrac(age, reference = "Younger")), 4L)
})

test_that("log_hr_from_ci() inverts hr_ci_from_log()", {
  x <- log_hr_from_ci(1.0, 0.8, 1.2, ci_level = 0.95)
  back <- hr_ci_from_log(x$yi, x$sei, ci_level = 0.95)
  round <- hr_ci_from_log(x$yi, x$sei, ci_level = 0.95)
  round2 <- log_hr_from_ci(round$hr, round$ci_lower, round$ci_upper, ci_level = 0.95)
  expect_equal(round2$yi, x$yi, tolerance = 1e-12)
  expect_equal(round2$sei, x$sei, tolerance = 1e-12)
})