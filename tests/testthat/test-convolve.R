test_that("convolve() recovers a normal sum", {
  skip_if_not_installed("bayesmeta")

  dens1 <- function(x, log = FALSE) stats::dnorm(x, mean = 1, sd = 2, log = log)
  dens2 <- function(x, log = FALSE) stats::dnorm(x, mean = 3, sd = 1, log = log)
  cdf1 <- function(x) stats::pnorm(x, mean = 1, sd = 2)
  cdf2 <- function(x) stats::pnorm(x, mean = 3, sd = 1)
  conv <- interact:::interactmeta_convolve_sum(dens1, dens2, cdf1, cdf2)
  # X+Y ~ N(4, sqrt(5))
  qs <- c(0.025, 0.5, 0.975)
  expect_equal(
    conv$quantile(qs),
    stats::qnorm(qs, mean = 4, sd = sqrt(5)),
    tolerance = 0.03
  )
  x <- 4
  expect_equal(
    conv$density(x),
    stats::dnorm(x, mean = 4, sd = sqrt(5)),
    tolerance = 0.01
  )
})

test_that("convolve() recovers a normal difference", {
  skip_if_not_installed("bayesmeta")

  dens_x <- function(x, log = FALSE) stats::dnorm(x, mean = 1, sd = 2, log = log)
  dens_y <- function(x, log = FALSE) stats::dnorm(x, mean = 3, sd = 1, log = log)
  cdf_x <- function(x) stats::pnorm(x, mean = 1, sd = 2)
  cdf_y <- function(x) stats::pnorm(x, mean = 3, sd = 1)
  dens_ny <- interact:::interactmeta_negate_density(dens_y)
  cdf_ny <- interact:::interactmeta_negate_cdf(cdf_y)
  conv <- interact:::interactmeta_convolve_sum(dens_x, dens_ny, cdf_x, cdf_ny)
  # X-Y ~ N(-2, sqrt(5))
  qs <- c(0.025, 0.5, 0.975)
  expect_equal(
    conv$quantile(qs),
    stats::qnorm(qs, mean = -2, sd = sqrt(5)),
    tolerance = 0.03
  )
})

test_that("scale density uses the absolute Jacobian", {
  dens <- function(x, log = FALSE) stats::dnorm(x, mean = 1, sd = 2, log = log)
  cdf <- function(x) stats::pnorm(x, mean = 1, sd = 2)

  dens_pos <- interact:::interactmeta_scale_density(dens, 2)
  z <- 3
  expect_equal(dens_pos(z), stats::dnorm(z / 2, mean = 1, sd = 2) / 2, tolerance = 1e-10)

  dens_neg <- interact:::interactmeta_scale_density(dens, -2)
  expect_equal(dens_neg(z), stats::dnorm(z / -2, mean = 1, sd = 2) / 2, tolerance = 1e-10)

  cdf_neg <- interact:::interactmeta_scale_cdf(cdf, -2)
  # Z = -2X, X ~ N(1, 4) => Z ~ N(-2, 16)
  expect_equal(cdf_neg(0), stats::pnorm(0, mean = -2, sd = 4), tolerance = 1e-8)
})

test_that("Gaussian mixture product weights are independent and normalized", {
  gm_x <- interact:::interactmeta_gm(c(0, 2), c(1, 1), c(0.25, 0.75))
  gm_y <- interact:::interactmeta_gm(c(-1, 1), c(0.5, 1.5), c(0.4, 0.6))
  gm <- interact:::interactmeta_gm_linear_sum(gm_x, 1, gm_y, -1)
  expect_equal(sum(gm$weight), 1)
  expect_equal(length(gm$weight), 4L)
  expect_equal(
    sort(gm$weight),
    sort(c(0.25, 0.75) %o% c(0.4, 0.6)),
    tolerance = 1e-12
  )
  expect_equal(
    interact:::interactmeta_gm_mean(gm),
    interact:::interactmeta_gm_mean(gm_x) - interact:::interactmeta_gm_mean(gm_y),
    tolerance = 1e-12
  )
  expect_equal(
    interact:::interactmeta_gm_var(gm),
    interact:::interactmeta_gm_var(gm_x) + interact:::interactmeta_gm_var(gm_y),
    tolerance = 1e-12
  )
})

test_that("mixture of two independent normals matches convolution", {
  skip_if_not_installed("bayesmeta")

  gm_x <- interact:::interactmeta_gm(1, 2, 1)
  gm_y <- interact:::interactmeta_gm(3, 1, 1)
  gm <- interact:::interactmeta_gm_linear_sum(gm_x, 1, gm_y, -1)
  sm <- interact:::interactmeta_gm_summary(gm)
  expect_equal(sm$mean, -2, tolerance = 1e-10)
  expect_equal(sm$sd, sqrt(5), tolerance = 1e-10)

  dens_x <- function(x, log = FALSE) stats::dnorm(x, 1, 2, log = log)
  dens_y <- function(x, log = FALSE) stats::dnorm(x, 3, 1, log = log)
  conv <- interact:::interactmeta_convolve_sum(
    dens_x,
    interact:::interactmeta_negate_density(dens_y),
    function(x) stats::pnorm(x, 1, 2),
    interact:::interactmeta_negate_cdf(function(x) stats::pnorm(x, 3, 1))
  )
  qs <- c(0.025, 0.5, 0.975)
  expect_equal(conv$quantile(qs), unname(sm$q[c(1, 2, 3)]), tolerance = 0.03)
})

test_that("exponentiated Gaussian-mixture summaries keep a finite sd", {
  mu <- -0.62
  s <- 0.15
  gm <- interact:::interactmeta_gm(mean = mu, sd = s, weight = 1)
  ln <- interact:::interactmeta_gm_lognormal_moments(gm)
  expect_equal(ln$mean, exp(mu + 0.5 * s^2), tolerance = 1e-12)
  expect_equal(
    ln$sd,
    exp(mu + 0.5 * s^2) * sqrt(exp(s^2) - 1),
    tolerance = 1e-12
  )

  draws <- data.frame(
    type = "Overall",
    effect = mu,
    .weight = 1,
    .cond_sd = s
  )
  fin <- interact:::interactmeta_summarize_effects(
    draws = draws,
    group_cols = "type",
    prob = 0.95,
    exponentiate = TRUE,
    fit = list(engine = "bayesmeta"),
    quantity = "combined",
    omega = 0.4
  )
  expect_true(is.finite(fin$summary$sd))
  expect_equal(fin$summary$mean, ln$mean, tolerance = 1e-8)
  expect_equal(fin$summary$sd, ln$sd, tolerance = 1e-8)
  expect_equal(fin$summary$estimate, exp(mu), tolerance = 1e-6)

  tab <- interact:::comparison_effects_print_df(fin$summary, digits = 2)
  expect_false(any(grepl("NA", unlist(tab), fixed = TRUE)))
})

test_that("CAMS reconstruction convolution matches analytic independent Gaussians", {
  skip_if_not_installed("bayesmeta")
  # Trial-level identities M - pi g and M + (1 - pi) g with independent
  # (M, g). Not population alpha + pi* theta.

  mu_M <- 0.2
  sd_M <- 0.5
  mu_g <- 0.4
  sd_g <- 0.3
  pi_star <- 0.35
  dens_M <- function(x, log = FALSE) stats::dnorm(x, mu_M, sd_M, log = log)
  dens_g <- function(x, log = FALSE) stats::dnorm(x, mu_g, sd_g, log = log)
  cdf_M <- function(x) stats::pnorm(x, mu_M, sd_M)
  cdf_g <- function(x) stats::pnorm(x, mu_g, sd_g)

  conv_A <- interact:::interactmeta_convolve_sum(
    dens_M,
    interact:::interactmeta_scale_density(dens_g, -pi_star),
    cdf_M,
    interact:::interactmeta_scale_cdf(cdf_g, -pi_star)
  )
  conv_B <- interact:::interactmeta_convolve_sum(
    dens_M,
    interact:::interactmeta_scale_density(dens_g, 1 - pi_star),
    cdf_M,
    interact:::interactmeta_scale_cdf(cdf_g, 1 - pi_star)
  )

  mu_A <- mu_M - pi_star * mu_g
  sd_A <- sqrt(sd_M^2 + (pi_star * sd_g)^2)
  mu_B <- mu_M + (1 - pi_star) * mu_g
  sd_B <- sqrt(sd_M^2 + ((1 - pi_star) * sd_g)^2)
  qs <- c(0.025, 0.5, 0.975)

  expect_equal(conv_A$quantile(qs), stats::qnorm(qs, mu_A, sd_A), tolerance = 0.03)
  expect_equal(conv_B$quantile(qs), stats::qnorm(qs, mu_B, sd_B), tolerance = 0.03)
  expect_equal(conv_A$density(mu_A), stats::dnorm(mu_A, mu_A, sd_A), tolerance = 0.02)
  expect_equal(conv_B$density(mu_B), stats::dnorm(mu_B, mu_B, sd_B), tolerance = 0.02)
  expect_equal(conv_A$cdf(mu_A), 0.5, tolerance = 0.02)
  expect_equal(conv_B$cdf(Inf), 1, tolerance = 0.01)
})

test_that("Monte Carlo reconstruction matches convolution quantiles", {
  skip_if_not_installed("bayesmeta")

  mu_M <- -0.3
  sd_M <- 0.4
  mu_g <- 0.25
  sd_g <- 0.2
  pi_star <- 0.6
  dens_M <- function(x, log = FALSE) stats::dnorm(x, mu_M, sd_M, log = log)
  dens_g <- function(x, log = FALSE) stats::dnorm(x, mu_g, sd_g, log = log)
  conv_A <- interact:::interactmeta_convolve_sum(
    dens_M,
    interact:::interactmeta_scale_density(dens_g, -pi_star),
    function(x) stats::pnorm(x, mu_M, sd_M),
    interact:::interactmeta_scale_cdf(
      function(x) stats::pnorm(x, mu_g, sd_g),
      -pi_star
    )
  )
  conv_B <- interact:::interactmeta_convolve_sum(
    dens_M,
    interact:::interactmeta_scale_density(dens_g, 1 - pi_star),
    function(x) stats::pnorm(x, mu_M, sd_M),
    interact:::interactmeta_scale_cdf(
      function(x) stats::pnorm(x, mu_g, sd_g),
      1 - pi_star
    )
  )

  set.seed(1)
  n <- 2e5
  M <- stats::rnorm(n, mu_M, sd_M)
  g <- stats::rnorm(n, mu_g, sd_g)
  muA <- M - pi_star * g
  muB <- M + (1 - pi_star) * g
  qs <- c(0.025, 0.5, 0.975)

  expect_equal(mean(muA), mu_M - pi_star * mu_g, tolerance = 0.01)
  expect_equal(stats::sd(muA), sqrt(sd_M^2 + (pi_star * sd_g)^2), tolerance = 0.01)
  expect_equal(as.numeric(stats::quantile(muA, qs)), conv_A$quantile(qs), tolerance = 0.02)
  expect_equal(as.numeric(stats::quantile(muB, qs)), conv_B$quantile(qs), tolerance = 0.02)
})

test_that("gm_linear_sum refuses two mixtures from the same CAMS block", {
  gm_a <- interact:::interactmeta_gm_tag_block(
    interact:::interactmeta_gm(0, 1, 1),
    "overall"
  )
  gm_th <- interact:::interactmeta_gm_tag_block(
    interact:::interactmeta_gm(1, 1, 1),
    "overall"
  )
  expect_error(
    interact:::interactmeta_gm_linear_sum(gm_a, 1, gm_th, 0.4),
    "independent mixtures"
  )
  gm_g <- interact:::interactmeta_gm_tag_block(
    interact:::interactmeta_gm(0.2, 0.5, 1),
    "interaction"
  )
  expect_silent(interact:::interactmeta_gm_linear_sum(gm_th, 1, gm_g, -1))
})

test_that("Var(alpha + pi theta) keeps Cov(alpha, theta); independent product does not", {
  w <- c(0.3, 0.7)
  a <- c(-0.2, 0.5)
  th <- c(-0.4, 0.8)
  va <- c(0.04, 0.09)
  vt <- c(0.16, 0.25)
  cat <- c(0.03, 0.08)
  pi <- 0.4

  mean_M <- a + pi * th
  var_M_node <- va + (pi^2) * vt + 2 * pi * cat
  gm_joint <- interact:::interactmeta_gm(mean_M, sqrt(pmax(var_M_node, 0)), w)
  var_joint <- interact:::interactmeta_gm_var(gm_joint)

  mu_a <- sum(w * a)
  mu_th <- sum(w * th)
  V_a <- sum(w * va) + sum(w * (a - mu_a)^2)
  V_th <- sum(w * vt) + sum(w * (th - mu_th)^2)
  Cov <- sum(w * cat) + sum(w * (a - mu_a) * (th - mu_th))
  expect_equal(var_joint, V_a + (pi^2) * V_th + 2 * pi * Cov, tolerance = 1e-12)
  expect_gt(abs(Cov), 1e-8)

  gm_a <- interact:::interactmeta_gm(a, sqrt(va), w)
  gm_th <- interact:::interactmeta_gm(th, sqrt(vt), w)
  gm_ind <- interact:::interactmeta_gm_linear_sum(gm_a, 1, gm_th, pi)
  expect_gt(
    abs(interact:::interactmeta_gm_var(gm_ind) - var_joint),
    1e-6
  )
})
