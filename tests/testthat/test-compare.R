test_that("comparison labels use inference when the analysis repeats", {
  bayes <- list(
    analysis = "cam",
    engine = "bayesmeta",
    method = "bayes",
    inference = "bayesian"
  )
  class(bayes) <- c("interactmeta", "cam")
  freq <- list(
    analysis = "cam",
    engine = "metafor",
    method = "REML",
    inference = "frequentist"
  )
  class(freq) <- c("interactmeta", "cam")
  da <- list(
    analysis = "da",
    engine = "metafor",
    method = "REML",
    inference = "frequentist"
  )
  class(da) <- "interactmeta_da"
  sw <- list(
    analysis = "swada",
    engine = "metafor",
    method = "REML",
    inference = "frequentist"
  )
  class(sw) <- "interactmeta_swada"

  expect_equal(
    interactmeta:::interactmeta_compare_labels(list(bayes, freq)),
    c("CAMS \u2014 bayesmeta", "CAMS \u2014 metafor (REML)")
  )
  expect_equal(
    interactmeta:::interactmeta_compare_labels(list(da, sw, bayes)),
    c("DA", "SWADA", "CAMS")
  )
  expect_equal(
    interactmeta:::interactmeta_inference_display(freq),
    "metafor (REML)"
  )
})

test_that("compare_* on two CAMS approaches labels and plots distinct approaches", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("ggplot2")

  fit_bayes <- suppressWarnings(cam(
    age,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))
  fit_freq <- cam(
    age,
    reference = "Younger",
    engine = "metafor",
    silent = 2
  )

  expect_silent(ix <- compare_interactions(fit_bayes, fit_freq, exponentiate = TRUE))
  expect_s3_class(ix, "interactmeta_compare_interactions")
  expect_equal(
    as.character(unique(ix$approach)),
    c("CAMS \u2014 bayesmeta", "CAMS \u2014 metafor (REML)")
  )
  expect_equal(nrow(ix), 2L)
  expect_true(all(is.finite(ix$estimate)))
  expect_true(all(is.finite(ix$ci.lb[ix$approach == "CAMS \u2014 bayesmeta"])))
  p_ix <- plot(ix)
  expect_s3_class(p_ix, "ggplot")

  expect_silent(sg <- compare_subgroups(fit_bayes, fit_freq, exponentiate = TRUE))
  expect_s3_class(sg, "interactmeta_compare_subgroups")
  expect_equal(nlevels(sg$approach), 2L)
  expect_equal(length(unique(as.character(sg$subgroup))), 2L)
  expect_s3_class(plot(sg), "ggplot")

  expect_silent(het <- compare_heterogeneity(fit_bayes, fit_freq))
  expect_s3_class(het, "interactmeta_compare_heterogeneity")
  expect_true("Overall" %in% as.character(het$parameter))
  expect_s3_class(plot(het), "ggplot")

  expect_silent(oi <- compare_ifrac(fit_bayes, fit_freq))
  expect_s3_class(oi, "interactmeta_compare_ifrac")
  expect_equal(nrow(oi), 2L)
  bayes_row <- oi[as.character(oi$approach) == "CAMS \u2014 bayesmeta", , drop = FALSE]
  freq_row <- oi[as.character(oi$approach) == "CAMS \u2014 metafor (REML)", , drop = FALSE]
  expect_true(is.finite(bayes_row$ci.lb))
  expect_true(is.na(freq_row$ci.lb))
  expect_s3_class(plot(oi), "ggplot")

  expect_silent(oi_g <- compare_ifrac(fit_bayes, fit_freq, ifrac = "interaction"))
  expect_s3_class(oi_g, "interactmeta_compare_ifrac")
  expect_equal(unique(oi_g$ifrac), "interaction")
  expect_equal(nrow(oi_g), 2L)
  expect_s3_class(plot(oi_g), "ggplot")

  expect_silent(comb <- compare_combined(fit_bayes, fit_freq, exponentiate = TRUE))
  expect_s3_class(comb, "interactmeta_compare_combined")
  expect_equal(nrow(comb), 2L)
  expect_s3_class(plot(comb), "ggplot")

  expect_error(
    compare_weights(fit_bayes, fit_freq),
    "compare_ifrac"
  )
})

test_that("mixed DA/SWADA/CAMS comparison objects remain plottable", {
  skip_if_not_installed("metafor")
  da_fit <- da(age, reference = "Younger")
  ad_fit <- ad(age, reference = "Younger")
  sw <- swada(age, reference = "Younger")
  cam <- cam(age, reference = "Younger", engine = "metafor", silent = 2)
  ix <- compare_interactions(da_fit, ad_fit, sw, cam)
  expect_s3_class(ix, "interactmeta_compare_interactions")
  expect_equal(as.character(ix$approach), c("DA", "AD", "SWADA", "CAMS"))
  expect_s3_class(plot(ix), "ggplot")
  w <- compare_weights(da_fit, ad_fit, sw)
  expect_s3_class(w, "interactmeta_compare_weights")
  expect_s3_class(plot(w), "ggplot")
})
