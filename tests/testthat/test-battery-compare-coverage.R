# compare_* tables, ggplot compare plots, and compare.R helpers (fast REML).

test_that("compare tables and plots cover main interactmeta_compare classes", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("ggplot2")
  b <- fixture_respire_constructions()
  reml <- b$cams
  ml <- fixture_reml_respire_cams("ML")

  het <- compare_heterogeneity(b$da, reml, component = "all")
  expect_s3_class(het, "interactmeta_compare_heterogeneity")
  p_het <- plot(het)
  expect_s3_class(p_het, "ggplot")

  w <- compare_weights(b$da, b$ad, b$swada)
  expect_true(all(is.finite(w$da_reference)))
  p_w <- plot(w)
  expect_s3_class(p_w, "ggplot")

  ifrac_int <- interactmeta:::compare_ifrac(reml, ml, ifrac = "interaction")
  expect_s3_class(ifrac_int, "interactmeta_compare_ifrac")
  p_if <- plot(ifrac_int)
  expect_s3_class(p_if, "ggplot")

  ifrac_ov <- interactmeta:::compare_ifrac(reml, ml, ifrac = "overall")
  expect_s3_class(ifrac_ov, "interactmeta_compare_ifrac")

  comb <- interactmeta:::compare_combined(reml, ml, type = "Overall", exponentiate = FALSE)
  expect_s3_class(comb, "interactmeta_compare_combined")
  p_cb <- plot(comb)
  expect_s3_class(p_cb, "ggplot")

  cam_cmp <- compare(reml, ml, type = "Overall", exponentiate = FALSE)
  expect_s3_class(cam_cmp, "interactmeta_compare_interactions")
  p_cam <- plot(cam_cmp)
  expect_s3_class(p_cam, "ggplot")
})

test_that("compare(sensitivity) modifier table is plottable", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("ggplot2")
  sens <- sensitivity(respire_age, reference = "<65", method = "REML", silent = 2)
  mod <- compare(sens, estimand = "modifiers", exponentiate = FALSE)
  expect_s3_class(mod, "interactmeta_compare_subgroups")
  p <- plot(mod)
  expect_s3_class(p, "ggplot")
})

test_that("plot_modifiers builds a ggplot from REML construction fits", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("ggplot2")
  b <- fixture_respire_constructions()
  p <- plot_modifiers(
    b$da,
    b$swada,
    b$cams,
    exponentiate = FALSE,
    include_observed = FALSE
  )
  expect_s3_class(p, "ggplot")
})