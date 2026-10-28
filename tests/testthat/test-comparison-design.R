test_that("prepare_ifrac_input() runs ifrac on escalc-like input", {
  es <- data.frame(
    study = rep("A", 2),
    level = c("Younger", "Older"),
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.02)
  )
  class(es) <- c("escalc", "data.frame")
  out <- interactmeta:::prepare_ifrac_input(
    es,
    reference = "Younger",
    trial = "trial",
    subgroup = "level",
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
    subgroup = "level",
    contrast = "all",
    contrast_sep = "-",
    big_var = 1e6,
    include_interaction = FALSE
  )
  expect_equal(out$ifrac_Older, 0.5)
})

test_that("resolve_comparison_design() builds formula for age-like ifrac output", {
  dat <- data.frame(
    study = rep("A", 2),
    level = c("Younger", "Older"),
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.02),
    reference = "Younger",
    ifrac_Older = 0.3,
    ifrac2_Older = c(-0.3, 0.7)
  )
  d <- interactmeta:::resolve_comparison_design(dat, reference = "Younger", trial_var = "study", subgroup_var = "level")
  expect_equal(d$ifrac_var, "ifrac_Older")
  expect_equal(d$ifrac2_var, "ifrac2_Older")
  expect_equal(d$effects, "ifrac_Older:level")
  expect_equal(d$n_nonref, 1L)
  expect_true(grepl("ifrac_1", d$formula))
  expect_true(grepl("\\|\\| study", d$formula))
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
    class = c("interactmeta", "list")
  )
  attr(fit, "interactmeta") <- list(
    ifrac_var = "ifrac_Older",
    ifrac2_var = "ifrac2_Older",
    trial_var = "study",
    subgroup_var = "level",
    effects = "ifrac_Older:subgroup"
  )
  st <- interactmeta:::interactmeta_settings(fit)
  expect_equal(st$ifrac_var, "ifrac_Older")
  expect_equal(st$effects, "ifrac_Older:subgroup")
})

test_that("comparison_model_formula() uses uncorrelated trial REs", {
  uncorr <- interactmeta:::comparison_model_formula(1L)
  expect_match(uncorr, "\\(1 \\+ ifrac2_1 \\|\\| study\\)")
})

test_that("default_comparison_prior() uses trial and ifrac2 coefficient names", {
  p <- interactmeta:::default_comparison_prior("ifrac2_Older", group = "trial")
  expect_equal(p$coef, c("Intercept", "ifrac2_Older"))
  expect_equal(p$group, c("trial", "trial"))
})

test_that("cams() returns interactmeta with restored names", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cams(
    dat,
    reference = "Younger",
    )
  expect_s3_class(fit, "interactmeta")
  expect_false(inherits(fit, "cam"))
  expect_equal(attr(fit, "interactmeta")$backend, "bayesmeta")
  expect_equal(as.character(fit$data$level), rep(c("Younger", "Older"), 2))
  expect_true("b_ifrac_Older" %in% names(fit$support))
  expect_true("b_subgroupOlder" %in% names(fit$support) ||
                "b_subgroupOlder" %in% gsub(" ", "", names(fit$support)))
  expect_equal(fit$data$yi, dat$yi)
})

test_that("cams() uses default seed for reproducible refits", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  args <- list(
    reference = "Younger",
        silent = 2
  )
  fit1 <- suppressWarnings(do.call(cams, c(list(dat), args)))
  fit2 <- suppressWarnings(do.call(cams, c(list(dat), args)))
  expect_equal(
    median(fit1$support$b_Intercept),
    median(fit2$support$b_Intercept)
  )

  fit_custom <- suppressWarnings(do.call(cams, c(list(dat), args, list(seed = 99L))))
  # DIRECT is deterministic; seed does not change the support.
  expect_equal(
    median(fit1$support$b_Intercept),
    median(fit_custom$support$b_Intercept)
  )
})

test_that("cams() rejects retired engine arguments", {
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  expect_error(
    cams(dat, reference = "Younger", use_precompiled = FALSE),
    "no longer supported"
  )
  expect_error(
    cams(dat, reference = "Younger", correlated_re = TRUE),
    "no longer supported"
  )
})
