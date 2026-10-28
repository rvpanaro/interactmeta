test_that("shrinkage(type = 'subgroups') matches modifiers(shrinkage = TRUE)", {
  fit <- fixture_cam_tiny_fit()
  me <- shrinkage(fit, type = "modifiers")
  alt <- modifiers(fit, shrinkage = TRUE)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_equal(me$summary, alt$summary)
})

test_that("shrinkage() type matches forestplot labels (no 'all')", {
  fit <- fixture_cam_tiny_fit()
  expect_s3_class(shrinkage(fit, type = "modifiers"), "interactmeta_subgroup_effects")
  expect_s3_class(shrinkage(fit, type = "interaction"), "interactmeta_interaction_effects")
  expect_s3_class(shrinkage(fit, type = "combined"), "interactmeta_combined_effect")
  expect_error(shrinkage(fit, type = "all"), "not supported")
})

test_that("shrinkage(component = ...) still works with a soft deprecation warning", {
  fit <- fixture_cam_tiny_fit()
  expect_warning(
    me <- shrinkage(fit, component = "subgroup"),
    "soft-deprecated"
  )
  expect_s3_class(me, "interactmeta_subgroup_effects")
})

test_that("modifiers(shrinkage = TRUE) returns one row per data row", {
  fit <- fixture_cam_tiny_fit()
  me <- modifiers(fit, shrinkage = TRUE)
  expect_true(me$shrinkage)
  expect_equal(unique(me$support$type), "Shrinkage")
  expect_equal(nrow(me$summary), 4L)
  expect_true(all(c("trial", "level") %in% names(me$summary)))
})

test_that("shrinkage(exponentiate = TRUE) runs on renamed CAMS fits", {
  fit <- fixture_cam_tiny_fit()
  me <- shrinkage(fit, exponentiate = TRUE)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_true(me$exponentiate)
  expect_true(all(me$support$effect > 0))
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
  fit <- cams(
    d,
        reference = "Lower",
    silent = 2
  )
  me <- shrinkage(fit)
  expect_s3_class(me, "interactmeta_subgroup_effects")
  expect_equal(nrow(me$summary), nrow(attr(fit, "interactmeta")$data_fit))
})

test_that("overall shrinkage is one draw vector per trial (like interaction)", {
  fit <- fixture_cam_tiny_fit()
  n_draw <- nrow(as.data.frame(fit))
  sg <- shrinkage(fit, type = "modifiers", summary = FALSE)
  int <- shrinkage(fit, type = "interaction", summary = FALSE)
  ov <- shrinkage(fit, type = "combined", summary = FALSE)

  expect_s3_class(ov, "interactmeta_combined_effect")
  expect_true(ov$shrinkage)
  expect_true("trial" %in% names(ov$support))
  expect_false("subgroup" %in% names(ov$support))
  expect_equal(nrow(ov$support), 2L * n_draw)
  expect_equal(as.integer(table(ov$support$trial)), c(n_draw, n_draw))
  expect_equal(nrow(int$support), nrow(ov$support))
  expect_equal(sort(unique(ov$support$trial)), c("A", "B"))

  # Legacy helper: mix at a fixed omega still works
  omega <- rep(0.35, n_draw)
  got <- interactmeta:::comparison_overall_shrinkage_draws(
    sg_draws = sg$support,
    omega = omega,
    reference = "Younger",
    data = fit$data,
    trial_var = "study",
    subgroup_var = "level"
  )
  for (tr in c("A", "B")) {
    ref <- sg$support$effect[sg$support$trial == tr & sg$support$level == "Younger"]
    non <- sg$support$effect[sg$support$trial == tr & sg$support$level == "Older"]
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
      ov$support$effect[ov$support$trial == tr],
      as.numeric(cmp$support$m[, tr]),
      tolerance = 1e-10
    )
  }
})

test_that("shrinkage attaches overall IF metadata when trials differ", {
  dat <- fixture_cam_tiny(
    ifrac_older = c(0.3, 0.3, 0.5, 0.5),
    ifrac2_older = c(-0.3, 0.7, -0.5, 0.5)
  )
  fit <- fixture_cam_tiny_fit(variant = "half_ifrac", dat = dat)
  weighted <- shrinkage(fit, type = "modifiers", summary = TRUE)
  expect_false(is.null(weighted$overall_ifrac))
  expect_equal(weighted$overall_ifrac_source, "fit")
  # Shrinkage is at study-specific π_j with native trial latent effects
  cmp <- shrinkage_components(fit)
  shr <- weighted$support$effect[
    weighted$support$trial == "A" & weighted$support$level == "Older"
  ]
  expect_equal(as.numeric(shr), as.numeric(cmp$support$y_B[, "A"]), tolerance = 1e-10)
})
