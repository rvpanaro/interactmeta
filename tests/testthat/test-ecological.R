test_that("ecological_slope = FALSE is not implemented for public engines", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  dat <- data.frame(
    trial = rep(c("RESPIRE1", "RESPIRE2"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.22),
    vi = c(0.04, 0.02, 0.05, 0.03),
    reference = "Younger"
  )
  expect_error(
    bayescam(dat, reference = "Younger", ecological_slope = FALSE),
    "not implemented"
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "metafor", ecological_slope = FALSE),
    "not implemented"
  )
})
