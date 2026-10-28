predict_cam_example_fit <- function() {
  dat <- data.frame(
    study = rep(paste0("T", 1:4), each = 2),
    level = factor(rep(c("Younger", "Older"), 4), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.28, 0.08, 0.22, 0.15, 0.30),
    vi = c(0.04, 0.02, 0.05, 0.03, 0.03, 0.04, 0.02, 0.025)
  )
  suppressWarnings(cams(
    dat,
    reference = "Younger",
        silent = 2
  ))
}

cams_if_draw_effects <- function(fit, pi) {
  post <- interactmeta:::comparison_support_df(fit)
  design <- interactmeta:::comparison_marginal_design(fit)
  sg <- interactmeta:::comparison_subgroup_coef_terms(fit, design)
  ref <- design$reference
  nonref <- setdiff(design$subgroup_levels, ref)[[1L]]
  beta_col <- grep("^b_ifrac", names(post), value = TRUE)
  beta_col <- beta_col[!grepl("ifrac2", beta_col)]
  gamma_col <- grep("^b_subgroup", names(post), value = TRUE)
  list(
    post = post,
    M = interactmeta:::comparison_marginal_overall_effect_draws(
      post, pi, design, fit, "Overall"
    ),
    mu_A = interactmeta:::comparison_marginal_fe_draws(
      post, pi, design, sg[[ref]]
    ),
    mu_B = interactmeta:::comparison_marginal_fe_draws(
      post, pi, design, sg[[nonref]]
    ),
    alpha = post$b_Intercept,
    b_ifrac = post[[beta_col]],
    gamma = post[[gamma_col]]
  )
}

test_that("predict_cams_by_if() uses central mixture intervals by default", {
  skip_if_not_installed("bayesmeta")

  fit <- predict_cam_example_fit()
  pi_star <- 0.35
  pred <- predict_cams_by_if(fit, if_grid = pi_star, level = 0.95)
  expect_s3_class(pred, "interactmeta_cam_if_pred")
  expect_equal(attr(pred, "interval.rule"), "mixture_equal_tail")
  expect_equal(sort(unique(pred)), c("Older", "Younger"))
  expect_true(all(c("ifrac", "estimate", "lower", "upper") %in% names(pred)))
  expect_true(all(c("overall_estimate", "overall_lower", "overall_upper") %in% names(pred)))

  d <- cams_if_draw_effects(fit, pi_star)
  w <- d$post$.weight
  sd_M <- interactmeta:::interactmeta_lincomb_cond_sd(d$post, 1, pi_star, 0)
  sd_A <- interactmeta:::interactmeta_lincomb_cond_sd(d$post, 1, pi_star, -pi_star)
  sd_B <- interactmeta:::interactmeta_lincomb_cond_sd(d$post, 1, pi_star, 1 - pi_star)
  qi_M <- interactmeta:::interactmeta_mixture_moments(
    mean = d$M, weight = w, cond_sd = sd_M, prob = 0.95
  )
  qi_A <- interactmeta:::interactmeta_mixture_moments(
    mean = d$mu_A, weight = w, cond_sd = sd_A, prob = 0.95
  )
  qi_B <- interactmeta:::interactmeta_mixture_moments(
    mean = d$mu_B, weight = w, cond_sd = sd_B, prob = 0.95
  )
  dA <- pred[pred == "Younger", ]
  dB <- pred[pred == "Older", ]
  expect_equal(dA$estimate, qi_A$estimate, tolerance = 1e-8)
  expect_equal(dA$lower, qi_A$lci, tolerance = 1e-8)
  expect_equal(dA$upper, qi_A$uci, tolerance = 1e-8)
  expect_equal(dB$estimate, qi_B$estimate, tolerance = 1e-8)
  expect_equal(dB$lower, qi_B$lci, tolerance = 1e-8)
  expect_equal(dB$upper, qi_B$uci, tolerance = 1e-8)
  expect_equal(dA$overall_estimate, qi_M$estimate, tolerance = 1e-8)
  expect_equal(dA$overall_lower, qi_M$lci, tolerance = 1e-8)
  expect_equal(dA$overall_upper, qi_M$uci, tolerance = 1e-8)
})

test_that("CAMS IF identities hold draw-by-draw, including pi = 0 and pi = 1", {
  skip_if_not_installed("bayesmeta")

  fit <- predict_cam_example_fit()
  for (pi in c(0, 0.25, 0.40, 0.55, 1)) {
    d <- cams_if_draw_effects(fit, pi)
    theta <- d$b_ifrac + d$gamma
    expect_equal(d$mu_B - d$mu_A, d$gamma, tolerance = 1e-12)
    expect_equal((1 - pi) * d$mu_A + pi * d$mu_B, d$M, tolerance = 1e-12)
    expect_equal(d$M, d$alpha + theta * pi, tolerance = 1e-12)
    expect_equal(d$M, d$alpha + (d$b_ifrac + d$gamma) * pi, tolerance = 1e-12)
    if (abs(pi) > 1e-12 && any(abs(d$gamma) > 1e-12)) {
      expect_false(isTRUE(all.equal(d$M, d$alpha + d$b_ifrac * pi, tolerance = 1e-8)))
    }
  }

  d0 <- cams_if_draw_effects(fit, 0)
  expect_equal(d0$mu_A, d0$alpha, tolerance = 1e-12)
  expect_equal(d0$mu_B, d0$alpha + d0$gamma, tolerance = 1e-12)
  expect_equal(d0$M, d0$alpha, tolerance = 1e-12)

  d1 <- cams_if_draw_effects(fit, 1)
  expect_equal(d1$mu_A, d1$alpha + d1$b_ifrac, tolerance = 1e-12)
  expect_equal(d1$mu_B, d1$alpha + d1$b_ifrac + d1$gamma, tolerance = 1e-12)
  expect_equal(d1$M, d1$alpha + d1$b_ifrac + d1$gamma, tolerance = 1e-12)
  expect_equal(d1$mu_B, d1$M, tolerance = 1e-12)
})

test_that("bayesmeta pointwise IF intervals are mixture quantiles of the draws", {
  skip_if_not_installed("bayesmeta")

  fit <- predict_cam_example_fit()
  pi <- 0.35
  d <- cams_if_draw_effects(fit, pi)
  post <- d$post
  expect_true(interactmeta:::interactmeta_has_cond_moments(post))
  expect_true(all(c(
    ".cond_var_alpha",
    ".cond_var_theta",
    ".cond_cov_alpha_theta",
    ".cond_sd_gamma"
  ) %in% names(post)))
  expect_true(any(post$.cond_var_alpha > 0 | post$.cond_var_theta > 0 |
    abs(post$.cond_cov_alpha_theta) > 0 | post$.cond_sd_gamma > 0))

  w <- post$.weight
  sd_M <- interactmeta:::interactmeta_lincomb_cond_sd(post, 1, pi, 0)
  sd_A <- interactmeta:::interactmeta_lincomb_cond_sd(post, 1, pi, -pi)
  sd_B <- interactmeta:::interactmeta_lincomb_cond_sd(post, 1, pi, 1 - pi)
  sd_M_formula <- sqrt(pmax(
    post$.cond_var_alpha + (pi^2) * post$.cond_var_theta +
      2 * pi * post$.cond_cov_alpha_theta,
    0
  ))
  expect_equal(sd_M, sd_M_formula, tolerance = 1e-12)
  expect_equal(
    sd_A^2,
    sd_M^2 + (pi^2) * (post$.cond_sd_gamma^2),
    tolerance = 1e-12
  )
  expect_equal(
    sd_B^2,
    sd_M^2 + ((1 - pi)^2) * (post$.cond_sd_gamma^2),
    tolerance = 1e-12
  )

  post_no_cov <- post
  post_no_cov$.cond_cov_alpha_theta <- 0
  sd_M_no_cov <- interactmeta:::interactmeta_lincomb_cond_sd(post_no_cov, 1, pi, 0)
  if (any(abs(post$.cond_cov_alpha_theta) > 1e-10)) {
    expect_gt(max(abs(sd_M^2 - sd_M_no_cov^2)), 1e-10)
  }

  qi_mix <- interactmeta:::comparison_draw_moments(
    d$M, prob = 0.95, weight = w, cond_sd = sd_M
  )
  qi_means <- interactmeta:::interactmeta_mixture_moments(
    mean = d$M, weight = w, cond_sd = rep(0, length(d$M)), prob = 0.95
  )
  # Conditional variance still separates the mixture from a quantile of means.
  expect_gt(
    abs(qi_mix$lci - qi_means$lci) + abs(qi_mix$uci - qi_means$uci),
    1e-8
  )

  pred <- predict_cams_by_if(fit, if_grid = pi, level = 0.95)
  row_A <- pred[pred == "Younger", ]
  expect_equal(row_A$overall_estimate, qi_mix$estimate, tolerance = 1e-8)
  expect_equal(row_A$overall_lower, qi_mix$lci, tolerance = 1e-8)
  expect_equal(row_A$overall_upper, qi_mix$uci, tolerance = 1e-8)
  qi_A <- interactmeta:::comparison_draw_moments(
    d$mu_A, prob = 0.95, weight = w, cond_sd = sd_A
  )
  expect_equal(row_A$lower, qi_A$lci, tolerance = 1e-8)
  expect_equal(row_A$upper, qi_A$uci, tolerance = 1e-8)
  expect_true(is.list(
    interactmeta:::interactmeta_bayesmeta_population_location(
      fit, pi = pi, gamma_scale = -pi, predictive = FALSE, prob = 0.95
    )
  ))
})

test_that("predict_cams_by_if(interval.type = central) uses equal-tail mixture rule", {
  skip_if_not_installed("bayesmeta")

  fit <- predict_cam_example_fit()
  pi_star <- 0.35
  pred <- predict_cams_by_if(
    fit,
    if_grid = pi_star,
    level = 0.95,
    interval.type = "central"
  )
  expect_equal(attr(pred, "interval.rule"), "mixture_equal_tail")
  expect_equal(attr(pred, "interval.type"), "central")

  d <- cams_if_draw_effects(fit, pi_star)
  w <- d$post$.weight
  sd_A <- interactmeta:::interactmeta_lincomb_cond_sd(d$post, 1, pi_star, -pi_star)
  qi_A <- interactmeta:::interactmeta_mixture_moments(
    mean = d$mu_A,
    weight = w,
    cond_sd = sd_A,
    prob = 0.95,
    interval.type = "central"
  )
  row_A <- pred[pred == "Younger", ]
  expect_equal(row_A$lower, qi_A$lci, tolerance = 1e-8)
  expect_equal(row_A$upper, qi_A$uci, tolerance = 1e-8)

  pred_short <- predict_cams_by_if(
    fit,
    if_grid = pi_star,
    level = 0.95,
    interval.type = "shortest"
  )
  expect_equal(attr(pred_short, "interval.rule"), "mixture_shortest")
})

test_that("plot_contribution() uses predict_cam_by_if()", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("ggplot2")

  fit <- predict_cam_example_fit()
  pred <- predict_cam_by_if(
    fit,
    if_grid = NULL,
    resolution = 11,
    level = 0.95,
    interval.type = "central"
  )

  p_prev <- plot_contribution(fit, resolution = 11, show_footnote = FALSE)
  built <- ggplot2::ggplot_build(p_prev)
  ribbon_i <- which(vapply(
    p_prev$layers,
    function(ly) inherits(ly$geom, "GeomRibbon"),
    logical(1)
  ))[[1L]]
  rib <- built$data[[ribbon_i]]
  expect_equal(
    sort(unique(rib$x)),
    sort(100 * unique(pred$ifrac)),
    tolerance = 1e-8
  )
})

test_that("plot_contribution() uses the IF grid and study-level points", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("ggplot2")

  fit <- predict_cam_example_fit()
  p <- plot_contribution(fit, resolution = 21, show_footnote = FALSE)
  expect_s3_class(p, "ggplot")
  geoms <- vapply(p$layers, function(ly) class(ly$geom)[[1L]], character(1))
  expect_true("GeomRibbon" %in% geoms)
  expect_true("GeomPoint" %in% geoms)
  expect_true("GeomLinerange" %in% geoms)

  built <- ggplot2::ggplot_build(p)
  point_layer <- which(vapply(
    p$layers,
    function(ly) inherits(ly$geom, "GeomPoint"),
    logical(1)
  ))
  expect_true(length(point_layer) >= 1L)
  xs <- unique(built$data[[point_layer[[1L]]]]$x)
  raw <- interactmeta:::resolve_bubbleplot_raw(
    fit,
    interactmeta:::resolve_bubbleplot_settings(fit)
  )
  if_frac_pct <- 100 * sort(unique(raw$ifrac_Older))
  expect_equal(sort(xs), if_frac_pct, tolerance = 1e-6)
  expect_gt(length(xs), 1L)
  pts <- built$data[[point_layer[[1L]]]]
  expect_true(all(is.na(pts$fill) | pts$fill %in% c("NA", "transparent")))
  expect_true(all(pts$shape == 21L))

  ribbon_layer <- which(vapply(
    p$layers,
    function(ly) inherits(ly$geom, "GeomRibbon"),
    logical(1)
  ))
  expect_true(length(ribbon_layer) >= 1L)
  n_x_ribbon <- length(unique(built$data[[ribbon_layer[[1L]]]]$x))
  expect_equal(n_x_ribbon, 21L)
})
