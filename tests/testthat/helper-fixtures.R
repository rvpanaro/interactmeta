# Shared test fixtures and tier helpers (sourced automatically by testthat).

slow <- function() {
  if (!identical(Sys.getenv("INTERACT_SLOW_TESTS"), "true")) {
    testthat::skip("Set INTERACT_SLOW_TESTS=true for slow plot/integration tests")
  }
}

.fixture_cache <- new.env(parent = emptyenv())

fixture_cam_tiny <- function(
    ifrac_older = c(0.3, 0.3, 0.55, 0.55),
    ifrac2_older = c(-0.3, 0.7, -0.55, 0.45),
    vi = 0.01) {
  data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(
      rep(c("Younger", "Older"), 2),
      levels = c("Younger", "Older")
    ),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = vi,
    reference = "Younger",
    ifrac_Older = ifrac_older,
    ifrac2_Older = ifrac2_older
  )
}

fixture_cam_tiny_interaction <- function() {
  fixture_cam_tiny(
    ifrac_older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_older = c(-0.3, 0.7, -0.4, 0.6)
  )
}

fixture_cam_tiny_fit <- function(method = "Bayesian",
                                 variant = "default",
                                 dat = NULL,
                                 engine = NULL,
                                 ...) {
  if (!is.null(engine)) {
    method <- if (identical(engine, "metafor")) "REML" else "Bayesian"
  }
  key <- paste("cam_tiny_v2", method, variant, sep = "_")
  if (exists(key, envir = .fixture_cache, inherits = FALSE)) {
    return(get(key, envir = .fixture_cache))
  }
  if (is.null(dat)) {
    dat <- if (identical(variant, "interaction")) {
      fixture_cam_tiny_interaction()
    } else {
      fixture_cam_tiny()
    }
  }
  if (interact:::interactmeta_cams_method_is_bayesian(method) ||
      identical(method, "bayes")) {
    testthat::skip_if_not_installed("bayesmeta")
  } else {
    testthat::skip_if_not_installed("metafor")
  }
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    method = method,
    silent = 2,
    precompute = TRUE,
    lazy_latent = FALSE,
    ...
  ))
  assign(key, fit, envir = .fixture_cache)
  fit
}

fixture_age_escalc <- function() {
  if (exists("age_escalc", envir = .fixture_cache, inherits = FALSE)) {
    return(get("age_escalc", envir = .fixture_cache))
  }
  testthat::skip_if_not_installed("metafor")
  es <- metafor::escalc(measure = "IRR", yi = age$yi, vi = age$vi, data = age)
  assign("age_escalc", es, envir = .fixture_cache)
  es
}

fixture_age_cam_fit <- function(...) {
  key <- "age_cam_bayesmeta_v2"
  if (exists(key, envir = .fixture_cache, inherits = FALSE)) {
    return(get(key, envir = .fixture_cache))
  }
  testthat::skip_if_not_installed("bayesmeta")
  testthat::skip_if_not_installed("metafor")
  es <- fixture_age_escalc()
  fit <- suppressWarnings(cams(
    es,
    reference = "Younger",
    silent = 2,
    precompute = TRUE,
    lazy_latent = FALSE,
    ...
  ))
  assign(key, fit, envir = .fixture_cache)
  fit
}

fixture_edss_escalc <- function() {
  if (exists("edss_escalc", envir = .fixture_cache, inherits = FALSE)) {
    return(get("edss_escalc", envir = .fixture_cache))
  }
  testthat::skip_if_not_installed("metafor")
  es <- metafor::escalc(measure = "IRR", yi = edss$yi, vi = edss$vi, data = edss)
  assign("edss_escalc", es, envir = .fixture_cache)
  es
}

fit_edss_engine <- function(engine, method = NULL, ...) {
  if (is.null(method)) {
    method <- if (identical(engine, "metafor")) "REML" else "Bayesian"
  }
  args <- list(
    x = fixture_edss_escalc(),
    reference = "Lower",
    method = method,
    silent = 2
  )
  extra <- list(...)
  suppressWarnings(do.call(cams, c(args, extra)))
}

fixture_cams_block <- function() {
  data.frame(
    study = rep(c("T1", "T2", "T3"), each = 2),
    level = factor(rep(c("Younger", "Older"), 3), levels = c("Younger", "Older")),
    yi = c(0.10, 0.30, 0.12, 0.28, 0.08, 0.24),
    vi = c(0.04, 0.04, 0.05, 0.03, 0.03, 0.05)
  )
}

cams_engine_example_data <- function() {
  data.frame(
    study = rep(paste0("T", 1:4), each = 2),
    level = factor(rep(c("Younger", "Older"), 4), levels = c("Younger", "Older")),
    yi = c(0.10, 0.25, 0.12, 0.28, 0.08, 0.22, 0.15, 0.30),
    vi = c(0.04, 0.02, 0.05, 0.03, 0.03, 0.04, 0.02, 0.025)
  )
}

cams_block_example_data <- function() {
  fixture_cams_block()
}

fixture_cams_native_shrinkage <- function(reverse_rows = FALSE,
                                           subgroup_levels = c("Younger", "Older")) {
  dat <- data.frame(
    study = rep(c("ADVANCE", "AFFIRM", "DEFINE", "TEMSO"), each = 2),
    level = factor(
      rep(c("Younger", "Older"), 4),
      levels = subgroup_levels
    ),
    yi = c(
      -0.3425, -0.6539,
      -1.2730, -0.7985,
      -0.9835, -0.4231,
      -0.4403, -0.3272
    ),
    vi = c(
      0.025, 0.068,
      0.018, 0.030,
      0.040, 0.047,
      0.040, 0.033
    ),
    stringsAsFactors = FALSE
  )
  if (isTRUE(reverse_rows)) {
    dat <- dat[rev(seq_len(nrow(dat))), , drop = FALSE]
  }
  dat
}

expect_origin_engine <- function(obj, engine) {
  org <- attr(obj, "interactmeta_origin", exact = TRUE)
  testthat::expect_false(is.null(org))
  testthat::expect_equal(org$fit_backend, engine)
}
