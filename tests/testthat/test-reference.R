test_that("interactmeta_reference is the single source of truth on fits and data", {
  dat <- data.frame(
    study = rep("A", 2),
    level = c("Younger", "Older"),
    yi = c(0.1, 0.2),
    vi = c(0.04, 0.01),
    reference = "Younger",
    ifrac_Older = 0.8,
    ifrac2_Older = c(-0.8, 0.2),
    stringsAsFactors = FALSE
  )
  expect_equal(interactmeta:::interactmeta_reference(dat), "Younger")
  expect_equal(interactmeta:::interactmeta_nonref(dat), "Older")
  expect_error(
    interactmeta:::interactmeta_reference(dat, reference = "Older"),
    "does not match the stored reference"
  )

  # Mismatched cams(reference=) vs existing design must error
  expect_error(
    interactmeta:::prepare_ifrac_input(dat, reference = "Older"),
    "does not match the stored reference"
  )
  out <- interactmeta:::prepare_ifrac_input(dat, reference = "Younger")
  expect_equal(unique(as.character(out$reference)), "Younger")
})

test_that("ifrac_<reference> design columns are rejected", {
  dat <- data.frame(
    trial = "A",
    subgroup = "Younger",
    yi = 0.1,
    vi = 0.01,
    reference = "Younger",
    ifrac_Younger = 0.7,
    ifrac2_Younger = -0.7,
    stringsAsFactors = FALSE
  )
  expect_error(
    interactmeta:::prepare_ifrac_input(dat, reference = "Younger"),
    "name the reference subgroup"
  )
})

test_that("resolve_forest_reference does not invent a reference", {
  sg <- data.frame(
    trial = "A",
    level = factor(c("Younger", "Older"), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2),
    vi = c(0.01, 0.02)
  )
  expect_error(
    interactmeta:::resolve_forest_reference(sg, sg),
    "Cannot resolve forest reference"
  )
  sg$reference <- "Younger"
  expect_equal(interactmeta:::resolve_forest_reference(sg, sg), "Younger")
})

test_that("fitted object stores and reuses the given reference", {
  skip_if_not_installed("bayesmeta")

  dat <- data.frame(
    study = rep(c("A", "B"), each = 2),
    level = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.4, 0.12, 0.28),
    vi = c(0.05, 0.05, 0.04, 0.03),
    n = 10
  )
  fit <- suppressWarnings(cams(
    dat,
        reference = "Younger"))
  expect_equal(interactmeta:::interactmeta_reference(fit), "Younger")
  expect_equal(attr(fit, "interactmeta")$reference, "Younger")
  expect_equal(interactmeta:::interactmeta_nonref(fit), "Older")

  cmp <- shrinkage_components(fit)
  expect_equal(cmp$reference, "Younger")
  expect_equal(cmp$nonref, "Older")
  # π is IF of Older (= nonref), reference weight is 1-π
  expect_equal(unname(cmp$pi[["A"]]), 0.5, tolerance = 1e-10)

  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = NULL,
    include_shrinkage = FALSE,
    include_overall = FALSE,
    include_prediction = FALSE,
    include_combined = FALSE
  )
  # Forest uses the same stored reference (no re-inference)
  expect_true("Younger" %in% panels$row_labels$colour_label |
    "Younger" %in% panels$row_labels$subheading |
    "Younger" %in% panels$row_labels$heading)
})
