cams_audit_example_data <- function() {
  data.frame(
    trial = rep(c("T1", "T2", "T3"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.10, 0.30, 0.12, 0.28, 0.08, 0.24),
    vi = c(0.04, 0.04, 0.05, 0.03, 0.03, 0.05)
  )
}

test_that("fit stores call, model_spec, priors, computation, assumptions, transformed_data", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_audit_example_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))

  expect_true(is.call(fit$call))
  expect_equal(as.character(fit$call[[1L]]), "cam")
  expect_s3_class(fit$model_spec, "interactmeta_model_spec")
  expect_equal(fit$model_spec$reference, "Younger")
  expect_equal(fit$model_spec$n_subgroups, 2L)
  expect_true(isTRUE(fit$model_spec$ecological_slope))
  expect_equal(fit$engine, "bayesmeta")
  expect_equal(fit$method, "bayes")
  expect_equal(fit$algorithm, "quadrature")
  expect_true(isTRUE(fit$computation$deterministic))
  expect_false(isTRUE(fit$computation$stochastic))
  expect_false(isTRUE(fit$computation$engine_args$convolve_derived))
  expect_equal(fit$assumptions$subgroup_covariance, "assumed_zero")
  expect_false(isTRUE(fit$assumptions$mu_includes_random_effects))
  expect_equal(fit$draws_semantics$row_type, "quadrature_node")
  expect_equal(fit$draws_semantics$weight_column, ".weight")
  expect_true(!is.null(fit$provenance$packages$interactmeta))
  expect_true(!is.null(fit$transformed_data$trial_blocks$overall))
  expect_true(!is.null(fit$priors$common))
  expect_true(!is.null(fit$priors$engine))
  expect_equal(fit$priors$tau$scale, 0.5)
})

test_that("audit() returns a structured object a third party can print", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_audit_example_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))
  a <- audit(fit)
  expect_s3_class(a, "interactmeta_audit")
  expect_equal(a$engine$engine, "bayesmeta")
  expect_equal(a$engine$algorithm, "quadrature")
  expect_true(!is.null(a$model$information_fraction))
  expect_true(!is.null(a$transformed_data$trial_blocks))
  expect_true(length(a$warnings) >= 1L)
  out <- paste(utils::capture.output(print(a)), collapse = "\n")
  expect_match(out, "CAMS statistical audit", fixed = TRUE)
  expect_match(out, "bayesmeta", fixed = TRUE)
  expect_match(out, "Half-Normal", fixed = TRUE)
  expect_match(out, "assumed_zero", fixed = TRUE)
  expect_match(
    out,
    "independent overall-block and interaction-block trial effects",
    fixed = TRUE
  )
  expect_false(grepl("||", out, fixed = TRUE))
  expect_false(grepl("ifrac2_", out, fixed = TRUE))
})

test_that("stored trial blocks match an independent hand calculation", {
  skip_if_not_installed("bayesmeta")
  dat <- cams_audit_example_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))
  blocks <- fit$transformed_data$trial_blocks
  t1 <- dat[dat$trial == "T1", , drop = FALSE]
  yA <- t1$yi[[1L]]
  yB <- t1$yi[[2L]]
  vA <- t1$vi[[1L]]
  vB <- t1$vi[[2L]]
  pi <- (1 / vB) / (1 / vA + 1 / vB)
  m <- (1 - pi) * yA + pi * yB
  g <- yB - yA
  vm <- 1 / (1 / vA + 1 / vB)
  vg <- vA + vB

  row <- blocks$overall[blocks$overall$trial == "T1", , drop = FALSE]
  int <- blocks$interaction[blocks$interaction$trial == "T1", , drop = FALSE]
  expect_equal(unname(row$yi), m, tolerance = 1e-12)
  expect_equal(unname(row$vi), vm, tolerance = 1e-12)
  expect_equal(unname(row$ifrac), pi, tolerance = 1e-12)
  expect_equal(unname(int$yi), g, tolerance = 1e-12)
  expect_equal(unname(int$vi), vg, tolerance = 1e-12)
  expect_equal(m - pi * g, yA, tolerance = 1e-12)
  expect_equal(m + (1 - pi) * g, yB, tolerance = 1e-12)
})

test_that("metafor provenance records frequentist method and unused priors", {
  skip_if_not_installed("metafor")
  dat <- cams_audit_example_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "metafor",
    method = "REML",
    silent = 2
  ))
  expect_equal(fit$engine, "metafor")
  expect_equal(fit$method, "REML")
  expect_equal(fit$algorithm, "optimization")
  expect_true(isTRUE(fit$computation$deterministic))
  expect_equal(fit$draws_semantics$row_type, "point_estimate")
  expect_false(isTRUE(fit$priors$engine$applied))
  a <- audit(fit)
  expect_match(
    paste(utils::capture.output(print(a)), collapse = "\n"),
    "REML",
    fixed = TRUE
  )
})
