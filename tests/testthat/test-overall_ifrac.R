test_that("overall_ifrac_compute() matches hand calculation for two trials", {
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = rep(c("Younger", "Older"), 2),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.6, 0.6)
  )
  tau <- c(0.1, 0.2)
  oi <- interact:::overall_ifrac_compute(
    data = dat,
    trial_var = "study",
    subgroup_var = "level",
    reference = "Younger",
    ifrac_var = "ifrac_Older",
    tau = tau
  )
  expect_length(oi, 2L)
  expect_true(all(is.finite(oi)))
  expect_true(oi[1] > 0.3 && oi[1] < 0.6)
})

test_that("ifrac() works on data frame with fixed tau", {
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = rep(c("Younger", "Older"), 2),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.6, 0.6)
  )
  oi <- ifrac(dat, tau = 0.15, reference = "Younger", ifrac = "overall")
  expect_s3_class(oi, "interactmeta_overall_ifrac")
  expect_type(as.numeric(oi), "double")
  expect_length(oi, 1L)
  summ <- ifrac(dat, tau = 0.15, reference = "Younger", ifrac = "overall", summary = TRUE)
  expect_equal(
    names(summ),
    c("estimate", "median", "mean", "sd", "lci", "uci", "prob")
  )
})

test_that("ifrac() works on interactmeta", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = rep(c("Younger", "Older"), 2),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cams(dat, reference = "Younger")
  class(fit) <- c("interactmeta", class(fit))
  oi <- ifrac(fit)
  expect_s3_class(oi, "interactmeta_overall_ifrac")
  expect_true(length(oi) > 1L)
})

test_that("interaction_ifrac_compute() uses inverse interaction variance weights", {
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = rep(c("Younger", "Older"), 2),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.01, 0.04, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.6, 0.6)
  )
  tau_gamma <- 0.1
  oi <- interact:::interaction_ifrac_compute(
    data = dat,
    trial_var = "study",
    subgroup_var = "level",
    reference = "Younger",
    ifrac_var = "ifrac_Older",
    tau_gamma = tau_gamma
  )
  den_a <- 0.01 + 0.04
  den_b <- 0.02 + 0.02
  w_a <- 1 / (tau_gamma^2 + den_a)
  w_b <- 1 / (tau_gamma^2 + den_b)
  expect_equal(oi, (0.3 * w_a + 0.6 * w_b) / (w_a + w_b))
})

test_that("overall-weighted IF uses overall-block Var(m)", {
  skip_if_not_installed("metafor")
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = rep(c("Younger", "Older"), 2),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.01, 0.04, 0.02, 0.02),
    stringsAsFactors = FALSE
  )
  fit <- cams(dat, reference = "Younger", method = "REML", silent = 2)
  st <- attr(fit, "interactmeta", exact = TRUE)
  tq <- cam_trial_quantities(fit)
  comp <- interact:::overall_ifrac_trial_components(
    data = st$data_fit,
    trial_var = st$trial_var,
    subgroup_var = st$subgroup_var,
    reference = st$reference,
    ifrac_var = st$ifrac_vars_std[[1L]]
  )
  idx <- match(comp$trial, tq$trial)
  expect_equal(comp$den_overall, as.numeric(tq$v_m[idx]), tolerance = 1e-10)
  tau <- 0.15
  w_tau <- 1 / (tau^2 + as.numeric(tq$v_m[idx]))
  omega_tau <- sum(comp$pi * w_tau) / sum(w_tau)
  des <- ifrac(dat, reference = "Younger")
  expect_equal(
    as.numeric(ifrac(des, tau = tau, reference = "Younger", ifrac = "overall")),
    omega_tau,
    tolerance = 1e-10
  )
})

test_that("interaction_ifrac() differs from ifrac() on data", {
  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = rep(c("Younger", "Older"), 2),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = c(0.01, 0.04, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.6, 0.6)
  )
  oi <- ifrac(dat, tau = 0.15, reference = "Younger", ifrac = "overall")
  ii <- ifrac(dat, ifrac = "interaction", tau_gamma = 0.15, reference = "Younger")
  expect_s3_class(ii, "interactmeta_overall_ifrac")
  expect_equal(attr(ii, "ifrac"), "interaction")
  expect_false(isTRUE(all.equal(as.numeric(oi), as.numeric(ii))))
})
