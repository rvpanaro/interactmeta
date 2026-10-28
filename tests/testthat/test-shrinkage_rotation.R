test_that("CAMS shrinkage is study-specific at pi_j", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.22),
    vi = c(0.04, 0.02, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.60, 0.60, 0.40, 0.40),
    ifrac2_Older = c(-0.60, 0.40, -0.40, 0.60),
    n = 10
  )
  fit <- suppressWarnings(cams(
    dat,
        reference = "Younger"))

  sg <- shrinkage(fit, type = "modifiers", format = "list", summary = FALSE)
  int <- shrinkage(fit, type = "interaction", format = "list", summary = FALSE)
  ov <- shrinkage(fit, type = "combined", format = "list", summary = FALSE)
  cmp <- shrinkage_components(fit)
  expect_equal(sort(names(cmp$pi)), c("A", "B"))
  expect_equal(unname(cmp$pi[["A"]]), 0.6, tolerance = 1e-12)
  expect_equal(cmp$evaluation, "trial")
  expect_true(length(cmp$omega) >= 1L)

  for (tr in c("A", "B")) {
    pi_j <- unname(cmp$pi[[tr]])
    y_A <- cmp$support$y_A[, tr]
    y_B <- cmp$support$y_B[, tr]
    m <- cmp$support$m[, tr]
    g <- cmp$support$g[, tr]

    # Rebuild at study-specific π_j
    expect_equal(y_A, m - pi_j * g, tolerance = 1e-10)
    expect_equal(y_B, m + (1 - pi_j) * g, tolerance = 1e-10)
    expect_equal(m, (1 - pi_j) * y_A + pi_j * y_B, tolerance = 1e-10)
    expect_equal(g, y_B - y_A, tolerance = 1e-10)

    # Public shrinkage() outputs match components
    expect_equal(
      ov$support$effect[ov$support$trial == tr],
      as.numeric(m),
      tolerance = 1e-10
    )
    expect_equal(
      int$support$effect[int$support$trial == tr],
      as.numeric(g),
      tolerance = 1e-10
    )
    expect_equal(
      sg$support$effect[sg$support$trial == tr & sg$support$level == "Younger"],
      as.numeric(y_A),
      tolerance = 1e-10
    )
  }
})

test_that("standardized differs from shrinkage when pi_j != omega", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.22),
    vi = c(0.04, 0.02, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.60, 0.60, 0.40, 0.40),
    ifrac2_Older = c(-0.60, 0.40, -0.40, 0.60),
    n = 10
  )
  fit <- suppressWarnings(cams(
    dat, reference = "Younger",
    ))
  shr <- shrinkage(fit, type = "modifiers", summary = FALSE, at = 0.5)
  std <- standardized_trial_effects(
    fit, type = "modifiers", summary = FALSE, at = 0.5
  )
  # Trial A has π=0.6 ≠ 0.5
  ya_shr <- shr$support$effect[
    shr$support$trial == "A" & shr$support$level == "Younger"
  ]
  ya_std <- std$support$effect[
    std$support$trial == "A" & std$support$level == "Younger"
  ]
  expect_gt(max(abs(ya_shr - ya_std)), 1e-6)
})

test_that("shrinkage and modifiers(shrinkage=TRUE) share draws", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.3, 0.12, 0.22),
    vi = c(0.04, 0.04, 0.05, 0.03),
    reference = "Younger",
    ifrac_Older = c(0.5, 0.5, 0.4, 0.4),
    ifrac2_Older = c(-0.5, 0.5, -0.4, 0.6),
    n = 10
  )
  fit <- suppressWarnings(cams(
    dat, reference = "Younger",
    ))
  oi <- ifrac(fit, ifrac = "overall")
  a <- shrinkage(fit, type = "modifiers", summary = TRUE, at = oi)
  b <- modifiers(fit, shrinkage = TRUE, summary = TRUE, at = oi)
  expect_equal(a$summary$estimate, b$summary$estimate, tolerance = 1e-10)
})
