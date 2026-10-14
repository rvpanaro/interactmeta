test_that("shrinkage(type = 'subgroups') matches subgroups(shrinkage = TRUE)", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- shrinkage(fit, type = "subgroups")
  alt <- subgroups(fit, shrinkage = TRUE)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_equal(me$summary, alt$summary)
})

test_that("shrinkage() type matches forestplot labels (no 'all')", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")

  expect_s3_class(shrinkage(fit, type = "subgroups"), "interactmeta_subgroup_effects")
  expect_s3_class(shrinkage(fit, type = "interaction"), "interactmeta_interaction_effects")
  expect_s3_class(shrinkage(fit, type = "combined"), "interactmeta_combined_effect")
  expect_error(shrinkage(fit, type = "all"), "not supported")
})

test_that("shrinkage(component = ...) still works with a soft deprecation warning", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  expect_warning(
    me <- shrinkage(fit, component = "subgroup"),
    "soft-deprecated"
  )
  expect_s3_class(me, "interactmeta_subgroup_effects")
})

test_that("subgroups(shrinkage = TRUE) returns one row per data row", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- subgroups(fit, shrinkage = TRUE)
  expect_true(me$shrinkage)
  expect_equal(unique(me$draws$type), "Shrinkage")
  expect_equal(nrow(me$summary), 4L)
  expect_true(all(c("trial", "subgroup") %in% names(me$summary)))
})

test_that("shrinkage(exponentiate = TRUE) runs on renamed CAMS fits", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- shrinkage(fit, exponentiate = TRUE)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_true(me$exponentiate)
  expect_true(all(me$draws$effect > 0))
})

test_that("shrinkage() works when ifrac_data includes derived rows", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")

  dat <- edss
  e <- metafor::escalc(measure = "IRR", yi = dat$yi, vi = dat$vi, data = dat)
  d <- ifrac(
    e,
    reference = "Lower",
    include_combined = TRUE,
    include_interaction = TRUE
  )
  fit <- cam(
    d,
    engine = "bayesmeta",
    reference = "Lower",
    silent = 2
  )
  me <- shrinkage(fit)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_equal(nrow(me$summary), nrow(attr(fit, "interactmeta")$data_fit))
})

test_that("overall shrinkage is one draw vector per trial (like interaction)", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_Older = c(-0.3, 0.7, -0.55, 0.45)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  n_draw <- nrow(as.data.frame(fit))
  sg <- shrinkage(fit, type = "subgroups", summary = FALSE)
  int <- shrinkage(fit, type = "interaction", summary = FALSE)
  ov <- shrinkage(fit, type = "combined", summary = FALSE)

  expect_s3_class(ov, "interactmeta_combined_effect")
  expect_true(ov$shrinkage)
  expect_true("trial" %in% names(ov$draws))
  expect_false("subgroup" %in% names(ov$draws))
  expect_equal(nrow(ov$draws), 2L * n_draw)
  expect_equal(as.integer(table(ov$draws$trial)), c(n_draw, n_draw))
  expect_equal(nrow(int$draws), nrow(ov$draws))
  expect_equal(sort(unique(ov$draws$trial)), c("A", "B"))

  # Legacy helper: mix at a fixed omega still works
  omega <- rep(0.35, n_draw)
  got <- interact:::comparison_overall_shrinkage_draws(
    sg_draws = sg$draws,
    omega = omega,
    reference = "Younger",
    data = fit$data,
    trial_var = "trial",
    subgroup_var = "subgroup"
  )
  for (tr in c("A", "B")) {
    ref <- sg$draws$effect[sg$draws$trial == tr & sg$draws$subgroup == "Younger"]
    non <- sg$draws$effect[sg$draws$trial == tr & sg$draws$subgroup == "Older"]
    expect_equal(
      got$effect[got$trial == tr],
      (1 - omega) * ref + omega * non,
      tolerance = 1e-10
    )
  }

  # Production overall uses the same π_j-based m_j as components
  cmp <- shrinkage_components(ov)
  for (tr in c("A", "B")) {
    expect_equal(
      ov$draws$effect[ov$draws$trial == tr],
      as.numeric(cmp$draws$m[, tr]),
      tolerance = 1e-10
    )
  }
})

test_that("shrinkage attaches overall IF metadata when trials differ", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.5, 0.5),
    ifrac2_Older = c(-0.3, 0.7, -0.5, 0.5)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  weighted <- shrinkage(fit, type = "subgroups", summary = TRUE)
  expect_false(is.null(weighted$overall_ifrac))
  expect_equal(weighted$overall_ifrac_source, "fit")
  # Shrinkage is at study-specific π_j with stored REs
  cmp <- shrinkage_components(fit)
  shr <- weighted$draws$effect[
    weighted$draws$trial == "A" & weighted$draws$subgroup == "Older"
  ]
  expect_equal(as.numeric(shr), as.numeric(cmp$draws$y_B[, "A"]), tolerance = 1e-10)
})
