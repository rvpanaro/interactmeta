test_that("prepare_ifrac_input() runs ifrac on escalc-like input", {
  es <- data.frame(
    trial = rep("A", 2),
    subgroup = c("Younger", "Older"),
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.02)
  )
  class(es) <- c("escalc", "data.frame")
  out <- interactmeta:::prepare_ifrac_input(
    es,
    reference = "Younger",
    trial = "trial",
    subgroup = "subgroup",
    contrast = "all",
    contrast_sep = "-",
    big_var = 1e6,
    include_interaction = FALSE
  )
  expect_true("ifrac_Older" %in% names(out))
  expect_true("ifrac2_Older" %in% names(out))
})

test_that("prepare_ifrac_input() passes through ifrac output", {
  dat <- data.frame(
    trial = "A",
    subgroup = "Older",
    yi = 0.1,
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = 0.5,
    ifrac2_Older = 0.5
  )
  out <- interactmeta:::prepare_ifrac_input(
    dat,
    reference = "Younger",
    trial = "trial",
    subgroup = "subgroup",
    contrast = "all",
    contrast_sep = "-",
    big_var = 1e6,
    include_interaction = FALSE
  )
  expect_equal(out$ifrac_Older, 0.5)
})

test_that("resolve_comparison_design() builds formula for age-like ifrac output", {
  dat <- data.frame(
    trial = rep("A", 2),
    subgroup = c("Younger", "Older"),
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.02),
    reference = "Younger",
    ifrac_Older = 0.3,
    ifrac2_Older = c(-0.3, 0.7)
  )
  d <- interactmeta:::resolve_comparison_design(dat, reference = "Younger", trial_var = "trial", subgroup_var = "subgroup")
  expect_equal(d$ifrac_var, "ifrac_Older")
  expect_equal(d$ifrac2_var, "ifrac2_Older")
  expect_equal(d$effects, "ifrac_Older:subgroup")
  expect_equal(d$n_nonref, 1L)
  expect_true(grepl("ifrac_1", d$formula))
  expect_true(grepl("\\|\\| trial", d$formula))
})

test_that("cam(engine = 'cmdstanr') is rejected", {
  es <- data.frame(
    trial = "A",
    subgroup = "Older",
    yi = 0.1,
    vi = 0.01
  )
  expect_error(
    cam(es, reference = "Younger", engine = "cmdstanr"),
    "cmdstanr"
  )
})

test_that("interactmeta_settings() reads comparison attribute", {
  fit <- structure(
    list(data = data.frame(
      yi = 1,
      vi = 0.01,
      subgroup = "Older",
      trial = "A",
      ifrac_Older = 0.4
    )),
    class = c("interactmeta", "cam", "list")
  )
  attr(fit, "interactmeta") <- list(
    ifrac_var = "ifrac_Older",
    ifrac2_var = "ifrac2_Older",
    trial_var = "trial",
    subgroup_var = "subgroup",
    effects = "ifrac_Older:subgroup"
  )
  st <- interactmeta:::interactmeta_settings(fit)
  expect_equal(st$ifrac_var, "ifrac_Older")
  expect_equal(st$effects, "ifrac_Older:subgroup")
})

test_that("comparison_model_formula() uses uncorrelated trial REs", {
  uncorr <- interactmeta:::comparison_model_formula(1L)
  expect_match(uncorr, "\\(1 \\+ ifrac2_1 \\|\\| trial\\)")
})

test_that("default_comparison_prior() uses trial and ifrac2 coefficient names", {
  p <- interactmeta:::default_comparison_prior("ifrac2_Older", group = "trial")
  expect_equal(p$coef, c("Intercept", "ifrac2_Older"))
  expect_equal(p$group, c("trial", "trial"))
})

test_that("cam() returns interactmeta with restored names", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta")
  expect_s3_class(fit, "interactmeta")
  expect_s3_class(fit, "cam")
  expect_equal(attr(fit, "interactmeta")$engine, "bayesmeta")
  expect_equal(as.character(fit$data$subgroup), rep(c("Younger", "Older"), 2))
  expect_true("b_ifrac_Older" %in% names(fit$draws))
  expect_true("b_subgroupOlder" %in% names(fit$draws) ||
                "b_subgroupOlder" %in% gsub(" ", "", names(fit$draws)))
  expect_equal(fit$data$yi, dat$yi)
})

test_that("cam() uses default seed for reproducible refits", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  args <- list(
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  )
  fit1 <- suppressWarnings(do.call(cam, c(list(dat), args)))
  fit2 <- suppressWarnings(do.call(cam, c(list(dat), args)))
  expect_equal(
    median(fit1$draws$b_Intercept),
    median(fit2$draws$b_Intercept)
  )

  fit_custom <- suppressWarnings(do.call(cam, c(list(dat), args, list(seed = 99L))))
  # DIRECT is deterministic; seed does not change the support.
  expect_equal(
    median(fit1$draws$b_Intercept),
    median(fit_custom$draws$b_Intercept)
  )
})

test_that("cam() rejects retired engine arguments", {
  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  expect_error(
    cam(dat, reference = "Younger", use_precompiled = FALSE),
    "no longer supported"
  )
  expect_error(
    cam(dat, reference = "Younger", correlated_re = TRUE),
    "no longer supported"
  )
})

test_that("detect_ifrac_var() finds single ifrac_* column", {
  dat <- data.frame(ifrac_Older = 1, ifrac2_Older = 0, yi = 1)
  expect_equal(interactmeta:::detect_ifrac_var(dat), "ifrac_Older")
})
