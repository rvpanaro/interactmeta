test_that("bayesmeta native combined and interactions match DIRECT grid", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  fit <- fixture_age_cam_fit()

  ns <- asNamespace("interact")
  unlock <- function(name) {
    if (bindingIsLocked(name, ns)) {
      unlockBinding(name, ns)
    }
  }
  lock_back <- function(name, fun) {
    assign(name, fun, envir = ns)
    lockBinding(name, ns)
  }
  old_comb <- get("interactmeta_try_native_combined", envir = ns)
  old_int <- get("interactmeta_try_native_interactions", envir = ns)
  on.exit({
    unlock("interactmeta_try_native_combined")
    unlock("interactmeta_try_native_interactions")
    lock_back("interactmeta_try_native_combined", old_comb)
    lock_back("interactmeta_try_native_interactions", old_int)
  }, add = TRUE)

  unlock("interactmeta_try_native_combined")
  unlock("interactmeta_try_native_interactions")
  assign("interactmeta_try_native_combined", function(...) NULL, envir = ns)
  assign("interactmeta_try_native_interactions", function(...) NULL, envir = ns)
  lockBinding("interactmeta_try_native_combined", ns)
  lockBinding("interactmeta_try_native_interactions", ns)
  grid_comb <- combined(fit, type = "Overall", summary = TRUE)
  grid_int <- interactions(fit, type = "Overall", summary = TRUE)

  unlock("interactmeta_try_native_combined")
  unlock("interactmeta_try_native_interactions")
  lock_back("interactmeta_try_native_combined", old_comb)
  lock_back("interactmeta_try_native_interactions", old_int)
  native_comb <- combined(fit, type = "Overall", summary = TRUE)
  native_int <- interactions(fit, type = "Overall", summary = TRUE)

  cols <- c("estimate", "mean", "sd", "lci", "uci")
  expect_equal(
    native_comb$summary[, cols],
    grid_comb$summary[, cols],
    tolerance = 1e-5
  )
  expect_equal(
    native_int$summary[, cols],
    grid_int$summary[, cols],
    tolerance = 1e-5
  )
})

test_that("Prediction still uses grid path for combined", {
  skip_if_not_installed("bayesmeta")
  dat <- data.frame(
    study = rep(c("A", "B", "C"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.18),
    vi = c(0.04, 0.05, 0.03, 0.04, 0.05, 0.06)
  )
  fit <- cams(dat, reference = "Younger", silent = 2)
  pred <- combined(fit, type = "Prediction", summary = TRUE)
  expect_gt(nrow(pred$support), 0L)
  expect_true(all(is.finite(pred$summary$lci)))
})
