# Fast frequentist CAMS paths (shared accessors without bayesmeta dependency).

test_that("REML tiny CAMS: effect accessors match Bayesian smoke shape", {
  skip_if_not_installed("metafor")
  fit <- fixture_cam_tiny_fit(method = "REML")

  sg <- modifiers(fit)
  expect_s3_class(sg, "interactmeta_subgroup_effects")
  expect_equal(nrow(sg$summary), 4L)

  int <- interactions(fit)
  expect_s3_class(int, "interactmeta_interaction_effects")
  expect_true(all(is.finite(int$summary$estimate)))

  cb <- combined(fit, type = "Overall", summary = TRUE)
  expect_s3_class(cb, "interactmeta_combined_effect")

  het <- heterogeneity(fit, summary = TRUE)
  expect_s3_class(het, "interactmeta_heterogeneity")
  expect_gt(nrow(het$summary), 0L)

  sh <- shrinkage(fit, type = "modifiers", summary = TRUE)
  expect_gt(nrow(sh$summary), 0L)
})

test_that("ML and REML tiny CAMS both print method label", {
  skip_if_not_installed("metafor")
  for (meth in c("REML", "ML")) {
    fit <- fixture_cam_tiny_fit(method = meth)
    out <- paste(capture.output(print(fit)), collapse = "\n")
    expect_match(out, meth, fixed = TRUE)
  }
})

test_that("heterogeneity() components on REML respire_age CAMS", {
  skip_if_not_installed("metafor")
  fit <- fixture_reml_respire_cams("REML")
  for (comp in c("all", "intercept", "slope")) {
    het <- heterogeneity(fit, component = comp, summary = TRUE)
    expect_s3_class(het, "interactmeta_heterogeneity")
    expect_gt(nrow(het$summary), 0L)
  }
  td <- tidy(fit)
  expect_s3_class(td, "tbl_df")
  expect_gt(nrow(td), 0L)
})

test_that("subgroup_long reshapes wide trial by subgroup columns", {
  wide <- data.frame(
    study = c("T1", "T2"),
    yi_Younger = c(0.1, 0.15),
    vi_Younger = c(0.04, 0.05),
    yi_Older = c(0.2, 0.25),
    vi_Older = c(0.03, 0.04),
    stringsAsFactors = FALSE
  )
  long <- subgroup_long(
    wide,
    subgroups = c("Younger", "Older"),
    id_cols = "study",
    values = c("yi", "vi")
  )
  expect_equal(nrow(long), 4L)
  expect_true(all(c("subgroup", "yi", "vi", "study") %in% names(long)))
  expect_setequal(as.character(unique(long$subgroup)), c("Older", "Younger"))
})

test_that("subgroup_long returns long input unchanged", {
  dat <- data.frame(
    study = "A",
    subgroup = "Younger",
    yi = 0.1,
    vi = 0.01,
    stringsAsFactors = FALSE
  )
  out <- subgroup_long(dat)
  expect_equal(out, dat)
})

test_that("heterogeneity_footnote and HR helpers run on REML CAMS", {
  skip_if_not_installed("metafor")
  fit <- fixture_reml_respire_cams("REML")
  het <- heterogeneity(fit, summary = TRUE)
  note <- heterogeneity_footnote(het, digits = 2)
  expect_true(is.character(note))
  expect_gt(nchar(paste(note, collapse = "\n")), 10L)

  x <- log_hr_from_ci(1.0, 0.8, 1.2, ci_level = 0.95)
  back <- hr_ci_from_log(x$yi, x$sei, ci_level = 0.95)
  expect_equal(back$hr, 1.0, tolerance = 1e-12)
  expect_true(all(is.finite(c(back$hr, back$ci_lower, back$ci_upper))))
})