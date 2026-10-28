test_that("cams() default is Bayesian bayesmeta", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cams(dat, reference = "Younger", silent = 2)
  fit_b <- cams(dat, reference = "Younger", method = "Bayesian", silent = 2)
  expect_s3_class(fit, "interactmeta")
  expect_false(inherits(fit, "cam"))
  expect_false(inherits(fit, "cams"))
  expect_equal(interactmeta_fit_backend(fit), "bayesmeta")
  expect_equal(fit$method, "Bayesian")
  expect_equal(interactmeta_fit_backend(fit_b), "bayesmeta")
  expect_equal(
    as.numeric(fit_b$support$.weight),
    as.numeric(fit$support$.weight),
    tolerance = 1e-10
  )
  expect_true(is.data.frame(fit$support))
})

test_that("combined() returns only the combined-effect class", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cams(dat, reference = "Younger", silent = 2)
  cb <- combined(fit, type = "Overall", summary = TRUE)
  expect_s3_class(cb, "interactmeta_combined_effect")
  expect_false(inherits(cb, "interactmeta_overall_effect"))
  expect_true(all(cb$summary$type %in% "Overall"))
  expect_true(is.data.frame(cb$support))
})

test_that("removed overall aliases are rejected", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.05, 0.03, 0.04)
  )
  fit <- cams(dat, reference = "Younger", silent = 2)
  expect_false(exists("overall", envir = asNamespace("interact"), inherits = FALSE))
  expect_false(exists("cam", envir = asNamespace("interact"), inherits = FALSE))
  expect_false(exists("bayescams", envir = asNamespace("interact"), inherits = FALSE))
  expect_false(exists("subgroups", envir = asNamespace("interact"), inherits = FALSE))
  expect_error(observed_effects(fit, type = "overall"), "arg")
  expect_equal(observed_effects(fit, type = "combined")$type, "combined")
})

test_that("method = REML and ML use metafor", {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  by_reml <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    method = "REML",
    silent = 2
  ))
  by_freq <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    method = "frequentist",
    silent = 2
  ))
  by_ml <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    method = "ML",
    silent = 2
  ))
  expect_equal(interactmeta_fit_backend(by_reml), "metafor")
  expect_equal(by_reml$method, "REML")
  expect_equal(interactmeta_fit_backend(by_freq), "metafor")
  expect_equal(by_freq$method, "REML")
  expect_equal(interactmeta_fit_backend(by_ml), "metafor")
  expect_equal(by_ml$method, "ML")
  expect_equal(nrow(by_reml$support), 1L)
})

test_that("Bayesian method aliases use bayesmeta", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  by_bayes <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    method = "bayes",
    silent = 2
  ))
  by_default <- cams(dat, reference = "Younger", silent = 2)
  expect_equal(interactmeta_fit_backend(by_bayes), "bayesmeta")
  expect_equal(by_bayes$method, "Bayesian")
  expect_equal(
    as.numeric(by_bayes$support$.weight),
    as.numeric(by_default$support$.weight),
    tolerance = 1e-10
  )
})

test_that("removed approach and engine arguments error", {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  expect_error(
    cams(dat, reference = "Younger", method = "REML", silent = 2),
    "approach"
  )
  expect_error(
    cams(dat, reference = "Younger", engine = "metafor", silent = 2),
    "engine"
  )
  expect_error(
    cams(dat, reference = "Younger", engine_args = list(), silent = 2),
    "engine_args"
  )
})

test_that("modifier names the column", {
  expect_equal(interactmeta_modifier_heading("dose"), "Dose")
  expect_equal(interactmeta_modifier_heading("Dose"), "Dose")
  expect_equal(interactmeta_modifier_heading("subgroup"), "Subgroup")
  expect_equal(interactmeta_resolve_modifier_column("dose", "subgroup"), "dose")
  expect_equal(interactmeta_resolve_modifier_column(NULL, NULL, "subgroup"), "subgroup")

  skip_if_not_installed("metafor")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    dose = factor(rep(c("Low", "High"), 3), levels = c("Low", "High")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cams(
    dat,
    reference = "Low",
    modifier = "dose",
    method = "REML",
    silent = 2
  )
  expect_equal(attr(fit, "interactmeta")$modifier, "dose")
  expect_equal(attr(fit, "interactmeta")$subgroup_var, "level")
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Dose levels:", fixed = TRUE)
  expect_match(out, "Low", fixed = TRUE)
  expect_match(out, "High", fixed = TRUE)
  expect_false(grepl("subgroup effects", out, fixed = TRUE))
  s <- summary(fit)
  sout <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(sout, "Dose-specific effects", fixed = TRUE)
  m_sum <- modifiers(fit, summary = TRUE)$summary
  expect_true(nrow(m_sum) >= 1L)
})
