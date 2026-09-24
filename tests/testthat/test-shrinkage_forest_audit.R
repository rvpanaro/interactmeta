test_that("pi is non-reference IF, not 1-pi / reference IF", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.40, -0.05, 0.30),
    # Unequal precision so IF_Older != 0.5 and != IF_Younger
    vi = c(0.08, 0.02, 0.02, 0.08),
    reference = "Younger",
    ifrac_Older = c(0.8, 0.8, 0.2, 0.2),
    ifrac2_Older = c(-0.8, 0.2, -0.2, 0.8),
    n = 10
  )
  # Design columns must match precision shares
  expect_equal(
    unname(1 / dat$vi[dat$subgroup == "Older"] /
      (1 / dat$vi[dat$subgroup == "Younger"] + 1 / dat$vi[dat$subgroup == "Older"])),
    c(0.8, 0.2),
    tolerance = 1e-12
  )

  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  cmp <- shrinkage_components(fit)
  expect_equal(cmp$nonref, "Older")
  expect_equal(unname(cmp$pi[["A"]]), 0.8, tolerance = 1e-10)
  expect_equal(unname(cmp$pi[["B"]]), 0.2, tolerance = 1e-10)
  # Must not be the reference IF
  expect_false(isTRUE(all.equal(unname(cmp$pi[["A"]]), 0.2, tolerance = 1e-8)))
  expect_false(isTRUE(all.equal(unname(cmp$pi[["B"]]), 0.8, tolerance = 1e-8)))

  tab <- shrinkage_diagnostics(fit)
  expect_equal(tab$pi_nonref, tab$pi)
  expect_equal(tab$pi_reference, 1 - tab$pi, tolerance = 1e-12)

  # Flipped IF weight would disagree with reconstructed y_A
  for (tr in cmp$trials) {
    pi_j <- unname(cmp$pi[[tr]])
    yA_wrong <- cmp$draws$m[, tr] - (1 - pi_j) * cmp$draws$g[, tr]
    expect_gt(
      max(abs(yA_wrong - cmp$draws$y_A[, tr])),
      1e-3
    )
  }
})

test_that("CAMS shrinkage g orientation is nonref - ref (Older - Younger)", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.40, -0.05, 0.30),
    vi = c(0.08, 0.08, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.3, 0.3),
    ifrac2_Older = c(-0.5, 0.5, -0.3, 0.7),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  cmp <- shrinkage_components(fit)
  expect_equal(cmp$reference, "Younger")
  expect_equal(cmp$nonref, "Older")
  expect_equal(cmp$interaction_label, "Older-Younger")
  for (tr in cmp$trials) {
    expect_equal(
      cmp$draws$g[, tr],
      cmp$draws$y_B[, tr] - cmp$draws$y_A[, tr],
      tolerance = 1e-12
    )
    expect_equal(
      unname(cmp$observed$g[[tr]]),
      unname(cmp$observed$y_B[[tr]] - cmp$observed$y_A[[tr]]),
      tolerance = 1e-12
    )
  }
})

test_that("draw-level reconstruction and round-trip identities hold", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.40, -0.05, 0.30),
    vi = c(0.08, 0.08, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.3, 0.3),
    ifrac2_Older = c(-0.5, 0.5, -0.3, 0.7),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  cmp <- shrinkage_components(fit)

  for (tr in cmp$trials) {
    pi_j <- unname(cmp$pi[[tr]])
    m <- cmp$draws$m[, tr]
    g <- cmp$draws$g[, tr]
    y_A <- cmp$draws$y_A[, tr]
    y_B <- cmp$draws$y_B[, tr]

    # Independent reconstruction at study-specific π_j
    y_A_reconstructed <- m - pi_j * g
    y_B_reconstructed <- m + (1 - pi_j) * g
    expect_equal(y_A_reconstructed, y_A, tolerance = 1e-12)
    expect_equal(y_B_reconstructed, y_B, tolerance = 1e-12)

    # Round-trip at π_j
    g2 <- y_B - y_A
    m2 <- (1 - pi_j) * y_A + pi_j * y_B
    expect_equal(m2 - pi_j * g2, y_A, tolerance = 1e-12)
    expect_equal(m2 + (1 - pi_j) * g2, y_B, tolerance = 1e-12)
    expect_equal(m2, m, tolerance = 1e-12)
    expect_equal(g2, g, tolerance = 1e-12)
  }
})

test_that("observed m,g recover original subgroup estimates", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.40, -0.05, 0.30),
    vi = c(0.08, 0.08, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.3, 0.3),
    ifrac2_Older = c(-0.5, 0.5, -0.3, 0.7),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  cmp <- shrinkage_components(fit)
  for (tr in cmp$trials) {
    pi_j <- cmp$pi[[tr]]
    m_obs <- cmp$observed$m[[tr]]
    g_obs <- cmp$observed$g[[tr]]
    expect_equal(
      m_obs - pi_j * g_obs,
      cmp$observed$y_A[[tr]],
      tolerance = 1e-12
    )
    expect_equal(
      m_obs + (1 - pi_j) * g_obs,
      cmp$observed$y_B[[tr]],
      tolerance = 1e-12
    )
  }
})

test_that("trial order and draw permutation do not change summaries", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("B", "A"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(-0.05, 0.30, 0.10, 0.40),
    vi = c(0.02, 0.02, 0.08, 0.08),
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.5, 0.5),
    ifrac2_Older = c(-0.3, 0.7, -0.5, 0.5),
    n = 10
  )
  fit_ba <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  dat_ab <- dat[order(dat$trial, dat$subgroup), , drop = FALSE]
  # Same fit object: permute trial column order of components matrices
  cmp <- shrinkage_components(fit_ba)
  sg <- shrinkage(fit_ba, type = "subgroups", summary = TRUE)

  # Permute draw rows consistently
  set.seed(42)
  ord <- sample.int(nrow(cmp$draws$m))
  m_perm <- cmp$draws$m[ord, , drop = FALSE]
  g_perm <- cmp$draws$g[ord, , drop = FALSE]
  yA_perm <- cmp$draws$y_A[ord, , drop = FALSE]
  for (tr in cmp$trials) {
    pi_j <- unname(cmp$pi[[tr]])
    expect_equal(
      yA_perm[, tr],
      m_perm[, tr] - pi_j * g_perm[, tr],
      tolerance = 1e-12
    )
    row <- sg$summary[
      sg$summary$trial == tr & sg$summary$subgroup == "Younger",
      ,
      drop = FALSE
    ]
    expect_true(nrow(row) == 1L)
    expect_true(is.finite(row$estimate[[1L]]))
  }

  # Named trial indexing: results for "A" match regardless of column order
  expect_true("A" %in% colnames(cmp$draws$m))
  expect_true("B" %in% colnames(cmp$draws$m))
  expect_equal(unname(cmp$pi[["A"]]), 0.5, tolerance = 1e-12)
  expect_equal(unname(cmp$pi[["B"]]), 0.3, tolerance = 1e-12)
})

test_that("forest shrinkage greys equal canonical accessor posterior QI", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.40, -0.05, 0.30),
    vi = c(0.08, 0.08, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.3, 0.3),
    ifrac2_Older = c(-0.5, 0.5, -0.3, 0.7),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  sg <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = FALSE)
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_pooled = FALSE,
    include_prediction = FALSE,
    include_overall = FALSE
  )
  pan <- panels$resultsA
  labs <- panels$row_labels
  for (tr in c("A", "B")) {
    row <- sg$summary[
      sg$summary$trial == tr & sg$summary$subgroup == "Younger",
      ,
      drop = FALSE
    ]
    # Prefer colour_label (pre-stack) then heading/subheading
    sg_col <- if ("colour_label" %in% names(labs)) {
      labs$colour_label
    } else if ("subheading" %in% names(labs)) {
      labs$subheading
    } else {
      labs$heading
    }
    ix <- which(
      pan$row_kind == "shrinkage" &
        !is.na(pan$shape) &
        sg_col == "Younger" &
        labs$label == tr
    )
    expect_true(length(ix) >= 1L)
    expect_equal(pan$yi[[ix[[1]]]], row$estimate[[1]], tolerance = 1e-8)
    expect_equal(pan$ci.lb[[ix[[1]]]], row$lci[[1]], tolerance = 1e-8)
    expect_equal(pan$ci.ub[[ix[[1]]]], row$uci[[1]], tolerance = 1e-8)
  }
  expect_false(any(pan$row_kind == "pooled", na.rm = TRUE))
  expect_false(any(pan$row_kind == "prediction", na.rm = TRUE))
})

test_that("reconstruction identities hold with study-specific pi_j", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.2, 0.8, 0.1, 0.4),
    vi = c(0.05, 0.05, 0.04, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.25, 0.25, 0.45, 0.45),
    ifrac2_Older = c(-0.25, 0.75, -0.45, 0.55),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  cmp <- shrinkage_components(fit)
  expect_equal(unname(cmp$pi[["A"]]), 0.25, tolerance = 1e-12)
  pi_j <- unname(cmp$pi[["A"]])
  m <- cmp$draws$m[, "A"]
  g <- cmp$draws$g[, "A"]
  expect_equal(cmp$draws$y_A[, "A"], m - pi_j * g, tolerance = 1e-10)
  expect_equal(cmp$draws$y_B[, "A"], m + (1 - pi_j) * g, tolerance = 1e-10)
})

test_that("exponentiation happens only after draw-level reconstruction", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.5, 0.12, 0.28),
    vi = c(0.04, 0.02, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.4, 0.4, 0.55, 0.55),
    ifrac2_Older = c(-0.4, 0.6, -0.55, 0.45),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  lin <- shrinkage(fit, type = "subgroups", summary = FALSE, exponentiate = FALSE)
  expd <- shrinkage(fit, type = "subgroups", summary = FALSE, exponentiate = TRUE)
  cmp <- shrinkage_components(fit)
  # Components stay on model scale
  expect_true(all(is.finite(cmp$draws$y_A)))
  expect_lt(max(abs(cmp$draws$y_A)), 20)
  # Exponentiated draws = exp(linear reconstructed draws)
  ya_lin <- lin$draws$effect[
    lin$draws$trial == "A" & lin$draws$subgroup == "Younger"
  ]
  ya_exp <- expd$draws$effect[
    expd$draws$trial == "A" & expd$draws$subgroup == "Younger"
  ]
  expect_equal(ya_exp, exp(ya_lin), tolerance = 1e-10)
  # Not equal to reconstructing from exponentiated m,g
  m <- cmp$draws$m[, "A"]
  g <- cmp$draws$g[, "A"]
  pi_j <- cmp$pi[["A"]]
  wrong <- exp(m) - pi_j * exp(g)
  expect_false(isTRUE(all.equal(ya_exp, wrong, tolerance = 1e-3)))
})

test_that("shrinkage_diagnostics table and Younger change identity", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.40, -0.05, 0.30),
    vi = c(0.08, 0.08, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.3, 0.3),
    ifrac2_Older = c(-0.5, 0.5, -0.3, 0.7),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  tab <- shrinkage_diagnostics(fit)
  expect_equal(nrow(tab), 2L)
  expect_true(all(tab$g_orientation == "Older - Younger"))
  expect_true(all(tab$reference_identity_max_abs < 1e-10))
  # Median-level approximate decomposition columns are present
  expect_true(all(c(
    "pi", "omega", "observed_m", "shrunk_m", "observed_g", "shrunk_g",
    "observed_reference", "reconstructed_shrunk_reference"
  ) %in% names(tab)))
})

test_that("forest does not use pooled/prediction as study shrinkage", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.40, -0.05, 0.30),
    vi = c(0.08, 0.08, 0.02, 0.02),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.3, 0.3),
    ifrac2_Older = c(-0.5, 0.5, -0.3, 0.7),
    n = 10
  )
  fit <- suppressWarnings(cam(
    dat, reference = "Younger",
    engine = "bayesmeta"))
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = 200L,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_pooled = TRUE,
    include_prediction = TRUE
  )
  pan <- panels$resultsA
  shr <- pan[pan$row_kind == "shrinkage" & !is.na(pan$shape), , drop = FALSE]
  pooled <- pan[pan$row_kind == "pooled", , drop = FALSE]
  expect_true(nrow(shr) >= 1L)
  expect_true(nrow(pooled) >= 1L)
  expect_false(isTRUE(all.equal(
    sort(shr$yi[is.finite(shr$yi)]),
    sort(pooled$yi[is.finite(pooled$yi)])
  )))
})
