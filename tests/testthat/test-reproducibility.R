test_that("hand-calculated same-IF MAP predictive example matches", {

  source <- "RESPIRE1"
  target <- "RESPIRE2"
  tau_m_draws <- c(0.1, 0.2)
  tau_g_draws <- c(0.05, 0.15)
  tau_m <- tau_m_draws
  tau_g <- tau_g_draws

  out <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0.2, se_m_source = 0.1,
    m_target = 0.5, se_m_target = 0.2,
    g_source = 0.05, se_g_source = 0.15,
    g_target = 0.25, se_g_target = 0.15,
    tau_m = tau_m, tau_g = tau_g,
    source = source, target = target
  )

  delta_m <- 0.5 - 0.2
  q_m <- delta_m^2 / (0.1^2 + 0.2^2 + 2 * tau_m_draws^2)
  p_m <- mean(stats::pchisq(q_m, df = 1, lower.tail = FALSE))

  delta_g <- 0.25 - 0.05
  q_g <- delta_g^2 / (0.15^2 + 0.15^2 + 2 * tau_g_draws^2)
  p_g <- mean(stats::pchisq(q_g, df = 1, lower.tail = FALSE))

  expect_equal(out$source, source)
  expect_equal(out$target, target)
  expect_equal(out$assumption, "same information fraction")
  expect_equal(out$method, "same_information_fraction")
  expect_equal(out$overall$observed_difference, delta_m)
  expect_null(out$overall$composition_adjustment)
  expect_null(out$overall$adjusted_difference)
  expect_equal(
    as.numeric(out$overall$q),
    q_m,
    tolerance = 1e-12
  )
  expect_equal(out$overall$p_value, p_m, tolerance = 1e-12)
  expect_equal(out$interaction$observed_difference, delta_g)
  expect_equal(
    as.numeric(out$interaction$q),
    q_g,
    tolerance = 1e-12
  )
  expect_equal(out$interaction$p_value, p_g, tolerance = 1e-12)
})

test_that("hand-calculated composition_decomposition diagnostic matches", {

  beta_draws <- c(0.4, 0.6)
  tau_m_draws <- c(0.1, 0.2)
  tau_g_draws <- c(0.05, 0.15)
  out <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0.2, se_m_source = 0.1,
    m_target = 0.5, se_m_target = 0.2,
    g_source = 0.05, se_g_source = 0.15,
    g_target = 0.25, se_g_target = 0.15,
    tau_m = tau_m_draws,
    tau_g = tau_g_draws,
    source = "RESPIRE1", target = "RESPIRE2",
    method = "composition_decomposition",
    pi_source = 0.25, pi_target = 0.55,
    beta = beta_draws
  )

  delta_m <- 0.5 - 0.2
  d_m <- delta_m - beta_draws * (0.55 - 0.25)
  q_m <- d_m^2 / (0.1^2 + 0.2^2 + 2 * tau_m_draws^2)
  p_m <- mean(stats::pchisq(q_m, df = 1, lower.tail = FALSE))
  expect_equal(out$assumption, "composition decomposition")
  expect_equal(
    as.numeric(out$overall$adjusted_difference),
    d_m,
    tolerance = 1e-12
  )
  expect_equal(out$overall$p_value, p_m, tolerance = 1e-12)
})

test_that("denominators include 2 * tau^2 (not one copy)", {

  se_m_s <- 0.1
  se_m_t <- 0.1
  se_g_s <- 0.2
  se_g_t <- 0.2
  tau <- 0.5
  out <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0, se_m_source = se_m_s,
    m_target = 1, se_m_target = se_m_t,
    g_source = 0, se_g_source = se_g_s,
    g_target = 1, se_g_target = se_g_t,
    tau_m = tau, tau_g = tau,
    source = "A", target = "B"
  )
  var_m_one_tau <- se_m_s^2 + se_m_t^2 + tau^2
  var_m_two_tau <- se_m_s^2 + se_m_t^2 + 2 * tau^2
  q_wrong <- 1^2 / var_m_one_tau
  q_right <- 1^2 / var_m_two_tau
  expect_equal(
    as.numeric(out$overall$q),
    as.numeric(q_right),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    as.numeric(out$overall$q),
    as.numeric(q_wrong),
    tolerance = 1e-8
  )))
  expect_equal(
    as.numeric(out$interaction$q),
    as.numeric(1^2 / (se_g_s^2 + se_g_t^2 + 2 * tau^2)),
    tolerance = 1e-12
  )
})

test_that("tau = 0 reduces denominators to sampling variances only", {

  se_m_s <- 0.2
  se_m_t <- 0.3
  se_g_s <- 0.1
  se_g_t <- 0.4
  tau0 <- 0
  out <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0.0, se_m_source = se_m_s,
    m_target = 0.5, se_m_target = se_m_t,
    g_source = 0.1, se_g_source = se_g_s,
    g_target = 0.4, se_g_target = se_g_t,
    tau_m = tau0, tau_g = tau0,
    source = "A", target = "B"
  )
  q_m_expected <- (0.5)^2 / (se_m_s^2 + se_m_t^2)
  q_g_expected <- (0.3)^2 / (se_g_s^2 + se_g_t^2)
  expect_equal(
    as.numeric(out$overall$q),
    q_m_expected,
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(out$interaction$q),
    q_g_expected,
    tolerance = 1e-12
  )
  expect_equal(out$overall$sampling_var_sum, se_m_s^2 + se_m_t^2)
  expect_equal(out$interaction$sampling_var_sum, se_g_s^2 + se_g_t^2)
})

test_that("primary overall test does not use beta or information fractions", {

  # Same m/g/se/tau; any beta/pi would change a composition-adjusted Q
  tau <- rep(0.1, 3)
  out <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0.1, se_m_source = 0.2,
    m_target = 0.7, se_m_target = 0.2,
    g_source = 0.2, se_g_source = 0.1,
    g_target = 0.5, se_g_target = 0.1,
    tau_m = tau, tau_g = tau,
    source = "A", target = "B"
  )
  expect_null(out$overall$composition_adjustment)
  expect_false("pi" %in% names(out))
  expect_false("beta" %in% names(out$overall))
  # Q uses raw delta_m only
  q_expected <- (0.6)^2 / (0.2^2 + 0.2^2 + 2 * 0.1^2)
  expect_equal(
    as.numeric(out$overall$q),
    rep(q_expected, 3),
    tolerance = 1e-12
  )
})

test_that("interaction calculation never uses pi or beta", {

  tau <- rep(0.1, 3)
  out <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0, se_m_source = 0.1,
    m_target = 0, se_m_target = 0.1,
    g_source = 0.2, se_g_source = 0.1,
    g_target = 0.5, se_g_target = 0.1,
    tau_m = tau, tau_g = tau,
    source = "A", target = "B"
  )
  expect_equal(out$interaction$observed_difference, 0.3)
  expect_false("composition_adjustment" %in% names(out$interaction))
  expect_false("pi" %in% names(out$interaction))
  expect_false("beta" %in% names(out$interaction))
})

test_that("reversing source/target flips D but not Q or p when tau is shared", {

  tau_m <- rep(0.2, 40)
  tau_g <- rep(0.1, 40)
  forward <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0.1, se_m_source = 0.2,
    m_target = 0.4, se_m_target = 0.25,
    g_source = 0.05, se_g_source = 0.15,
    g_target = 0.35, se_g_target = 0.18,
    tau_m = tau_m, tau_g = tau_g,
    source = "A", target = "B"
  )
  reverse <- interactmeta:::comparison_reproducibility_compute(
    m_source = 0.4, se_m_source = 0.25,
    m_target = 0.1, se_m_target = 0.2,
    g_source = 0.35, se_g_source = 0.18,
    g_target = 0.05, se_g_target = 0.15,
    tau_m = tau_m, tau_g = tau_g,
    source = "B", target = "A"
  )

  expect_equal(reverse$overall$observed_difference, -forward$overall$observed_difference)
  expect_equal(
    reverse$interaction$observed_difference,
    -forward$interaction$observed_difference
  )
  expect_equal(
    as.numeric(reverse$overall$q),
    as.numeric(forward$overall$q),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(reverse$interaction$q),
    as.numeric(forward$interaction$q),
    tolerance = 1e-12
  )
  expect_equal(reverse$overall$p_value, forward$overall$p_value, tolerance = 1e-12)
  expect_equal(reverse$interaction$p_value, forward$interaction$p_value, tolerance = 1e-12)
})

make_two_study_dat <- function() {
  data.frame(
    trial = rep(c("RESPIRE1", "RESPIRE2"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.22),
    vi = c(0.04, 0.02, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.30, 0.30, 0.55, 0.55),
    ifrac2_Older = c(-0.30, 0.70, -0.55, 0.45),
    n = 10
  )
}

fit_source_only <- function(dat, seed = 1L) {
  suppressWarnings(cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger",
    seed = seed
  ))
}

test_that("a fit containing both source and target studies is rejected", {
  skip_if_not_installed("bayesmeta")

  dat <- make_two_study_dat()
  fit_both <- fit_source_only(dat)
  target_dat <- dat[dat$trial == "RESPIRE2", , drop = FALSE]

  expect_error(
    reproducibility(
      fit_both,
      target_data = target_dat,
      source = "RESPIRE1",
      target = "RESPIRE2"
    ),
    "source study only"
  )
  expect_error(
    reproducibility(
      fit_both,
      target_data = target_dat,
      source = "RESPIRE1",
      target = "RESPIRE2"
    ),
    "must not contribute"
  )
})

test_that("one-trial source CAMS is not available on public engines", {
  skip_if_not_installed("bayesmeta")
  dat <- make_two_study_dat()
  expect_error(
    fit_source_only(dat[dat$trial == "RESPIRE1", , drop = FALSE])
  )
})
