freq_cams_data <- function() {
  data.frame(
    study = rep(paste0("T", 1:4), each = 2),
    level = factor(rep(c("Younger", "Older"), 4), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.28, 0.08, 0.22, 0.15, 0.30),
    vi = c(0.04, 0.02, 0.05, 0.03, 0.03, 0.04, 0.02, 0.025)
  )
}

fit_metafor_freq <- function() {
  skip_if_not_installed("metafor")
  cams(freq_cams_data(), reference = "Younger", method = "REML", silent = 2)
}

fit_bayesmeta_freq <- function() {
  skip_if_not_installed("bayesmeta")
  suppressWarnings(cams(
    freq_cams_data(),
    reference = "Younger",
        silent = 2
  ))
}

test_that("metafor print uses CI, not CrI or posterior median", {
  fit <- fit_metafor_freq()
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "REML model results", fixed = TRUE)
  expect_match(out, "Intercept", fixed = TRUE)
  expect_match(out, "Ecological_Older", fixed = TRUE)
  expect_match(out, "estimate", fixed = TRUE)
  expect_false(grepl("median", out, fixed = TRUE))
  expect_match(out, "confidence interval", ignore.case = TRUE)
  expect_false(grepl("CrI", out, fixed = TRUE))
  expect_false(grepl("credible", out, ignore.case = TRUE))
  expect_false(grepl("marginal posterior summary", out, fixed = TRUE))
  expect_false(grepl("posterior draws:", out, fixed = TRUE))

  int_out <- paste(capture.output(print(interactions(fit, type = "Overall"))), collapse = "\n")
  expect_false(grepl("Treatment-by-subgroup interaction", int_out, fixed = TRUE))
  expect_false(grepl("at overall information fraction", int_out, fixed = TRUE))
  expect_match(int_out, "estimate", fixed = TRUE)
  expect_match(int_out, "Overall 95% CI", fixed = TRUE)
  expect_false(grepl("confidence intervals\\)", int_out))
  expect_false(grepl("median", int_out))
  expect_false(grepl("CrI", int_out, fixed = TRUE))
  expect_false(grepl("credible", int_out, ignore.case = TRUE))

  sg_out <- paste(capture.output(print(modifiers(fit, type = "Overall"))), collapse = "\n")
  expect_match(sg_out, "Overall 95% CI", fixed = TRUE)
  expect_false(grepl("Interaction-weighted IF:", sg_out, fixed = TRUE))
  expect_false(grepl("CrI", sg_out, fixed = TRUE))
  expect_false(grepl("median", sg_out))

  comb_out <- paste(capture.output(print(combined(fit, type = "Overall"))), collapse = "\n")
  expect_match(comb_out, "Overall 95% CI", fixed = TRUE)
  expect_false(grepl("Reporting IF:", comb_out, fixed = TRUE))
  expect_false(grepl("CrI", comb_out, fixed = TRUE))

  oi_out <- paste(capture.output(print(ifrac(fit))), collapse = "\n")
  expect_match(oi_out, "IF", fixed = TRUE)
  expect_match(oi_out, "interaction-weighted", fixed = TRUE)
  expect_match(oi_out, "estimate", fixed = TRUE)
  expect_false(grepl("Overall 95% CI", oi_out, fixed = TRUE))
  expect_false(grepl("Interaction-weighted information fraction", oi_out, fixed = TRUE))
  expect_match(
    oi_out,
    paste0("at fitted interaction heterogeneity ", "\u03c4", "_", "\u03b3"),
    fixed = TRUE
  )
  expect_false(grepl("CrI", oi_out, fixed = TRUE))
  expect_false(grepl("credible", oi_out, ignore.case = TRUE))

  oi_ov <- paste(capture.output(print(ifrac(fit, ifrac = "overall"))), collapse = "\n")
  expect_match(oi_ov, "combined-weighted", fixed = TRUE)
  expect_match(
    oi_ov,
    paste0("at fitted heterogeneity ", "\u03c4"),
    fixed = TRUE
  )

  sum_out <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(sum_out, "Overall 95% CI", fixed = TRUE)
  expect_match(sum_out, "Prediction 95% CI", fixed = TRUE)
  expect_false(grepl("95% CrI", sum_out, fixed = TRUE))
  expect_match(
    sum_out,
    paste0("at fitted interaction heterogeneity ", "\u03c4", "_", "\u03b3"),
    fixed = TRUE
  )
  expect_equal(
    length(gregexpr(
      paste0("at fitted interaction heterogeneity ", "\u03c4", "_", "\u03b3"),
      sum_out,
      fixed = TRUE
    )[[1L]]),
    1L
  )
})

test_that("frequentist IF uses fitted tau / tau_gamma", {
  fit <- fit_metafor_freq()
  tau_m <- sqrt(max(as.numeric(fit$metafor$overall$tau2), 0))
  tau_g <- sqrt(max(as.numeric(fit$metafor$interaction$tau2), 0))
  st <- attr(fit, "interactmeta", exact = TRUE)
  d <- st$data_fit
  ifrac_var <- st$ifrac_vars_std[[1L]]

  oi <- ifrac(fit, ifrac = "overall")
  expect_length(as.numeric(oi), 1L)
  expect_equal(
    as.numeric(oi),
    interact:::overall_ifrac_compute(
      data = d,
      trial_var = st$trial_var,
      subgroup_var = st$subgroup_var,
      reference = st$reference,
      ifrac_var = ifrac_var,
      tau = tau_m
    ),
    tolerance = 1e-12
  )

  oi_g <- ifrac(fit, ifrac = "interaction")
  expect_length(as.numeric(oi_g), 1L)
  expect_equal(
    as.numeric(oi_g),
    interact:::interaction_ifrac_compute(
      data = d,
      trial_var = st$trial_var,
      subgroup_var = st$subgroup_var,
      reference = st$reference,
      ifrac_var = ifrac_var,
      tau_gamma = tau_g
    ),
    tolerance = 1e-12
  )
})

test_that("ifrac(fit, ifrac='overall') matches original subgroup labels after recoding", {
  skip_if_not_installed("metafor")
  toy <- data.frame(
    study = c(rep("Both1", 2), rep("Both2", 2), "Bonly"),
    level = c("A", "B", "A", "B", "B"),
    yi = c(-0.2, -0.05, -0.3, -0.1, -0.05),
    vi = c(0.04, 0.04, 0.03, 0.05, 0.01)
  )
  des <- ifrac(toy, reference = "A", contrast = NULL)
  fit <- cams(toy, reference = "A", method = "REML", silent = 2)
  tau_m <- sqrt(max(as.numeric(fit$metafor$overall$tau2), 0))
  from_des <- as.numeric(ifrac(
    des, tau = tau_m, reference = "A", ifrac_var = "ifrac_B", ifrac = "overall"
  ))
  from_fit <- as.numeric(ifrac(fit, ifrac = "overall"))
  expect_equal(from_fit, from_des, tolerance = 1e-10)
  expect_lt(from_fit, 0.9)
})

test_that("frequentist M(pi), mu_A, mu_B use the rma covariance including Cov(alpha, theta)", {
  fit <- fit_metafor_freq()
  pi <- 0.4
  ov <- fit$metafor$overall
  int <- fit$metafor$interaction
  alpha <- as.numeric(ov$beta[[1L]])
  theta <- as.numeric(ov$beta[[2L]])
  gamma <- as.numeric(int$beta[[1L]])
  vb <- as.matrix(ov$vb)
  vg <- as.numeric(int$vb[[1L]])

  m <- alpha + pi * theta
  mu_a <- m - pi * gamma
  mu_b <- m + (1 - pi) * gamma
  expect_equal(mu_b - mu_a, gamma, tolerance = 1e-12)
  expect_equal((1 - pi) * mu_a + pi * mu_b, m, tolerance = 1e-12)

  var_m <- as.numeric(t(c(1, pi)) %*% vb %*% c(1, pi))
  var_a <- var_m + (pi^2) * vg
  var_b <- var_m + ((1 - pi)^2) * vg
  z <- stats::qnorm(0.975)

  sg <- modifiers(fit, type = "Overall", at = pi)
  comb <- combined(fit, type = "Overall", at = pi)
  intn <- interactions(fit, type = "Overall")
  ref <- sg$summary$level == "Younger"
  non <- sg$summary$level == "Older"
  expect_equal(sg$summary$estimate[ref], mu_a, tolerance = 1e-8)
  expect_equal(sg$summary$estimate[non], mu_b, tolerance = 1e-8)
  expect_equal(comb$summary$estimate[[1L]], m, tolerance = 1e-8)
  expect_equal(intn$summary$estimate[intn$summary$type == "Overall"], gamma, tolerance = 1e-8)
  expect_equal(sg$summary$sd[ref], sqrt(var_a), tolerance = 1e-8)
  expect_equal(sg$summary$sd[non], sqrt(var_b), tolerance = 1e-8)
  expect_equal(comb$summary$sd[[1L]], sqrt(var_m), tolerance = 1e-8)
  expect_equal(sg$summary$lci[ref], mu_a - z * sqrt(var_a), tolerance = 1e-8)
  expect_equal(sg$summary$uci[ref], mu_a + z * sqrt(var_a), tolerance = 1e-8)
})

test_that("exponentiate keeps frequentist CI interpretation", {
  fit <- fit_metafor_freq()
  raw <- interactions(fit, type = "Overall", exponentiate = FALSE)
  expd <- interactions(fit, type = "Overall", exponentiate = TRUE)
  pooled <- raw$summary$type == "Overall"
  expect_equal(
    expd$summary$estimate[pooled],
    exp(raw$summary$estimate[pooled]),
    tolerance = 1e-10
  )
  expect_equal(
    expd$summary$lci[pooled],
    exp(raw$summary$lci[pooled]),
    tolerance = 1e-10
  )
  expect_equal(
    expd$summary$uci[pooled],
    exp(raw$summary$uci[pooled]),
    tolerance = 1e-10
  )
  out <- paste(capture.output(print(expd)), collapse = "\n")
  expect_match(out, "Overall 95% CI", fixed = TRUE)
  expect_false(grepl("CrI", out, fixed = TRUE))
})

test_that("bayesmeta print keeps CrI and posterior summaries", {
  fit <- fit_bayesmeta_freq()
  int_out <- paste(capture.output(print(interactions(fit, type = "Overall"))), collapse = "\n")
  expect_false(grepl("Treatment-by-subgroup interaction", int_out, fixed = TRUE))
  expect_match(int_out, "median", fixed = TRUE)
  expect_match(int_out, "Overall 95% CrI", fixed = TRUE)
  expect_false(grepl("95% CI", int_out, fixed = TRUE))

  fit_out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(fit_out, "credible", ignore.case = TRUE)
  expect_false(grepl("Wald confidence", fit_out, fixed = TRUE))

  oi <- ifrac(fit)
  expect_gt(length(as.numeric(oi)), 1L)
})

test_that("CAMS identities hold for both inferential approaches", {
  pi_star <- 0.35
  check_ids <- function(fit, exact) {
    sg <- modifiers(fit, type = "Overall", at = pi_star)
    comb <- combined(fit, type = "Overall", at = pi_star)
    intn <- interactions(fit, type = "Overall")
    sg_ov <- sg$summary[sg$summary$type == "Overall", , drop = FALSE]
    mu_a <- sg_ov$estimate[sg_ov == "Younger"][[1L]]
    mu_b <- sg_ov$estimate[sg_ov == "Older"][[1L]]
    m <- comb$summary$estimate[comb$summary$type == "Overall"][[1L]]
    g <- intn$summary$estimate[intn$summary$type == "Overall"][[1L]]
    tol <- if (isTRUE(exact)) 1e-10 else 0.05
    expect_equal(as.numeric(mu_b - mu_a), as.numeric(g), tolerance = tol)
    expect_equal(
      as.numeric((1 - pi_star) * mu_a + pi_star * mu_b),
      as.numeric(m),
      tolerance = tol
    )
  }
  check_ids(fit_metafor_freq(), exact = TRUE)
  check_ids(fit_bayesmeta_freq(), exact = FALSE)
})
