test_that("CAMS shrinkage is study-specific at pi_j with stored REs", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.22),
    vi = c(0.04, 0.02, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.60, 0.60, 0.40, 0.40),
    ifrac2_Older = c(-0.60, 0.40, -0.40, 0.60),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger"))

  sg <- shrinkage(fit, type = "subgroups", format = "list", summary = FALSE)
  int <- shrinkage(fit, type = "interaction", format = "list", summary = FALSE)
  ov <- shrinkage(fit, type = "combined", format = "list", summary = FALSE)
  cmp <- shrinkage_components(fit)
  expect_equal(sort(names(cmp$pi)), c("A", "B"))
  expect_equal(unname(cmp$pi[["A"]]), 0.6, tolerance = 1e-12)
  expect_equal(cmp$evaluation, "trial")
  expect_true(length(cmp$omega) >= 1L)

  prep <- interact:::comparison_prepare_marginal_draws(fit, ifrac = "overall")
  post <- prep$post
  alpha <- as.numeric(post$b_Intercept)
  beta <- as.numeric(post$b_ifrac_Older)
  gamma <- as.numeric(post$b_subgroupOlder)

  for (tr in c("A", "B")) {
    pi_j <- unname(cmp$pi[[tr]])
    y_A <- cmp$draws$y_A[, tr]
    y_B <- cmp$draws$y_B[, tr]
    m <- cmp$draws$m[, tr]
    g <- cmp$draws$g[, tr]

    u_col <- paste0("r_trial[", tr, ",Intercept]")
    v_col <- paste0("r_trial[", tr, ",ifrac2_Older]")
    has_re <- u_col %in% names(post)
    u <- if (has_re) as.numeric(post[[u_col]]) else rep(0, nrow(post))
    v <- if (v_col %in% names(post)) as.numeric(post[[v_col]]) else rep(0, nrow(post))

    if (has_re) {
      expect_equal(
        as.numeric(y_A),
        alpha + u + beta * pi_j - pi_j * v,
        tolerance = 1e-8
      )
      expect_equal(
        as.numeric(y_B),
        alpha + u + beta * pi_j + (1 - pi_j) * v + gamma,
        tolerance = 1e-8
      )
      expect_equal(as.numeric(g), gamma + v, tolerance = 1e-8)
    }

    # Rebuild at study-specific π_j
    expect_equal(y_A, m - pi_j * g, tolerance = 1e-10)
    expect_equal(y_B, m + (1 - pi_j) * g, tolerance = 1e-10)
    expect_equal(m, (1 - pi_j) * y_A + pi_j * y_B, tolerance = 1e-10)
    expect_equal(g, y_B - y_A, tolerance = 1e-10)

    if (has_re) {
      ya_lp <- interact:::comparison_effect_linpred_draws(
        post = post,
        evaluation_if = pi_j,
        design = prep$design,
        subgroup_label = "Younger",
        subgroup_coef = prep$sg_coef[["Younger"]],
        re = "stored",
        fit = fit,
        trial = tr
      )
      expect_equal(as.numeric(y_A), as.numeric(ya_lp), tolerance = 1e-8)
    }

    # Public shrinkage() outputs match components
    expect_equal(
      ov$draws$effect[ov$draws$trial == tr],
      as.numeric(m),
      tolerance = 1e-10
    )
    expect_equal(
      int$draws$effect[int$draws$trial == tr],
      as.numeric(g),
      tolerance = 1e-10
    )
    expect_equal(
      sg$draws$effect[sg$draws$trial == tr & sg$draws$subgroup == "Younger"],
      as.numeric(y_A),
      tolerance = 1e-10
    )
  }
})

test_that("standardized and study-specific coincide when pi_j equals omega", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.3, 0.12, 0.22),
    vi = c(0.04, 0.04, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.4, 0.4),
    ifrac2_Older = c(-0.5, 0.5, -0.4, 0.6),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  # Identity is for stored trial REs evaluated at π_j vs ω. Public engines use
  # native trial latents for shrinkage, which are not the FE linpred at ω.
  post <- as.data.frame(fit)
  skip_if_not(
    any(grepl("^r_trial\\[", names(post))),
    "coincide identity requires stored trial random effects"
  )
  # Force overall IF to the trial π so estimands coincide
  shr <- shrinkage(fit, type = "subgroups", summary = FALSE, contribution = 0.5)
  std <- standardized_trial_effects(
    fit, type = "subgroups", summary = FALSE, contribution = 0.5
  )
  # Trial A has π = 0.5 = omega, so study-specific and standardized coincide there.
  expect_equal(
    shr$draws$effect[shr$draws$trial == "A"],
    std$draws$effect[std$draws$trial == "A"],
    tolerance = 1e-10
  )

  cmp_shr <- shrinkage_components(fit, contribution = 0.5)
  cmp_std <- interact:::comparison_shrinkage_components_list(
    interact:::comparison_shrinkage_rotation_bundle(
      fit, overall_ifrac = 0.5, evaluation = "overall"
    )
  )
  expect_equal(cmp_shr$draws$y_A[, "A"], cmp_std$draws$y_A[, "A"], tolerance = 1e-10)
  expect_equal(cmp_shr$draws$y_B[, "A"], cmp_std$draws$y_B[, "A"], tolerance = 1e-10)
})

test_that("standardized differs from shrinkage when pi_j != omega", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.22),
    vi = c(0.04, 0.02, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.60, 0.60, 0.40, 0.40),
    ifrac2_Older = c(-0.60, 0.40, -0.40, 0.60),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  shr <- shrinkage(fit, type = "subgroups", summary = FALSE, contribution = 0.5)
  std <- standardized_trial_effects(
    fit, type = "subgroups", summary = FALSE, contribution = 0.5
  )
  # Trial A has π=0.6 ≠ 0.5
  ya_shr <- shr$draws$effect[
    shr$draws$trial == "A" & shr$draws$subgroup == "Younger"
  ]
  ya_std <- std$draws$effect[
    std$draws$trial == "A" & std$draws$subgroup == "Younger"
  ]
  expect_gt(max(abs(ya_shr - ya_std)), 1e-6)
})

test_that("shrinkage and subgroups(shrinkage=TRUE) share draws", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.3, 0.12, 0.22),
    vi = c(0.04, 0.04, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.4, 0.4),
    ifrac2_Older = c(-0.5, 0.5, -0.4, 0.6),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  oi <- ifrac(fit, ifrac = "overall")
  a <- shrinkage(fit, type = "subgroups", summary = TRUE, contribution = oi)
  b <- subgroups(fit, shrinkage = TRUE, summary = TRUE, contribution = oi)
  expect_equal(a$summary$estimate, b$summary$estimate, tolerance = 1e-10)
})
