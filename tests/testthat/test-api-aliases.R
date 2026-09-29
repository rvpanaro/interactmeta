test_that("cam() matches cams() on the same data", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  )
  fit_alias <- cams(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  )
  expect_s3_class(fit, "interactmeta")
  expect_s3_class(fit, "cam")
  expect_s3_class(fit, "cams")
  expect_equal(fit$engine, "bayesmeta")
  expect_equal(fit$method, "bayes")
  expect_equal(fit_alias$engine, fit$engine)
  expect_equal(
    as.numeric(fit_alias$draws$.weight),
    as.numeric(fit$draws$.weight),
    tolerance = 1e-10
  )
})

test_that("combined() is a numerical alias of overall()", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cam(dat, reference = "Younger", engine = "bayesmeta", silent = 2)
  ov <- overall(fit, type = "Pooled", summary = TRUE)
  cb <- combined(fit, type = "Pooled", summary = TRUE)
  expect_s3_class(cb, "interactmeta_combined_effect")
  expect_s3_class(cb, "interactmeta_overall_effect")
  expect_s3_class(ov, "interactmeta_combined_effect")
  expect_equal(ov$summary$estimate, cb$summary$estimate, tolerance = 1e-12)
  expect_equal(ov$summary$lci, cb$summary$lci, tolerance = 1e-12)
  expect_equal(ov$summary$uci, cb$summary$uci, tolerance = 1e-12)
  expect_equal(ov$draws$effect, cb$draws$effect, tolerance = 1e-12)
})

test_that("type = 'combined' matches type = 'overall' for observed effects", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.04, 0.05, 0.03, 0.04)
  )
  fit <- cam(dat, reference = "Younger", engine = "bayesmeta", silent = 2)
  a <- observed_effects(fit, type = "overall")
  b <- observed_effects(fit, type = "combined")
  expect_equal(a$estimate, b$estimate, tolerance = 1e-12)
  expect_equal(a$variance, b$variance, tolerance = 1e-12)
  expect_equal(a$type, "overall")
  expect_equal(b$type, "overall")
})

test_that("method = REML without engine uses metafor", {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  by_method <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    method = "REML",
    silent = 2
  ))
  by_approach <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    approach = "REML",
    silent = 2
  ))
  by_ml <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    method = "ML",
    silent = 2
  ))
  expect_equal(by_method$engine, "metafor")
  expect_equal(by_method$method, "REML")
  expect_equal(by_approach$engine, "metafor")
  expect_equal(by_approach$method, "REML")
  expect_equal(by_ml$engine, "metafor")
  expect_equal(by_ml$method, "ML")
  expect_error(
    cam(dat, reference = "Younger", engine = "bayesmeta", method = "REML"),
    "method = \"bayes\""
  )
})

test_that("approach = bayes uses bayesmeta when engine is omitted", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  by_approach <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    approach = "bayes",
    silent = 2
  ))
  by_method <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    method = "bayes",
    silent = 2
  ))
  by_alias <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayes",
    silent = 2
  ))
  expect_equal(by_approach$engine, "bayesmeta")
  expect_equal(by_approach$method, "bayes")
  expect_equal(by_method$engine, "bayesmeta")
  expect_equal(by_alias$engine, "bayesmeta")
})

test_that("cam() without method defaults to bayesmeta DIRECT", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- suppressWarnings(cam(dat, reference = "Younger", silent = 2))
  expect_equal(fit$engine, "bayesmeta")
  expect_equal(fit$method, "bayes")
})

test_that("cam(engine = 'cmdstanr') is rejected", {
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  expect_error(
    cam(dat, reference = "Younger", engine = "cmdstanr", silent = 2),
    "removed"
  )
})

test_that("modifier names the column and aliases the extractors", {
  expect_identical(modifiers, subgroups)
  expect_identical(plot_modifiers, plot_subgroups)
  expect_identical(compare_modifiers, compare_subgroups)
  expect_equal(interactmeta_modifier_heading("dose"), "Dose")
  expect_equal(interactmeta_modifier_heading("Dose"), "Dose")
  expect_equal(interactmeta_modifier_heading("subgroup"), "Subgroup")
  expect_equal(interactmeta_resolve_modifier_column("dose", "subgroup"), "dose")
  expect_equal(interactmeta_resolve_modifier_column(NULL, NULL, "subgroup"), "subgroup")

  skip_if_not_installed("metafor")
  dat <- data.frame(
    trial = rep(c("A", "B", "C"), each = 2),
    dose = factor(rep(c("Low", "High"), 3), levels = c("Low", "High")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cam(
    dat,
    reference = "Low",
    modifier = "dose",
    engine = "metafor",
    silent = 2
  )
  expect_equal(attr(fit, "interactmeta")$modifier, "dose")
  expect_equal(attr(fit, "interactmeta")$subgroup_var, "subgroup")
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Dose:", fixed = TRUE)
  expect_match(out, "Low", fixed = TRUE)
  expect_match(out, "High", fixed = TRUE)
  expect_false(grepl("subgroup effects", out, fixed = TRUE))
  s <- summary(fit)
  sout <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(sout, "Dose effects", fixed = TRUE)
})
