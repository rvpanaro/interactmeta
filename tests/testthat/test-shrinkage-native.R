# Engine-native CAMS shrinkage (bayesmeta DIRECT / metafor BLUP).
# Public engines expose shrinkage through their native latent-effect objects.

parse_ifrac_display <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  as.numeric(sub("%", "", x, fixed = TRUE)) / 100
}

forest_shrinkage_rows <- function(fit, group = "subgroup") {
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_overall = FALSE,
    include_prediction = FALSE,
    include_combined = FALSE,
    group = group
  )
  pan <- panels$resultsA
  labs <- panels$row_labels
  sg <- as.character(labs$colour_label)
  trial <- as.character(labs$label)
  ix <- which(pan$row_kind == "shrinkage" & !is.na(pan$shape))
  data.frame(
    trial = trial[ix],
    subgroup = sg[ix],
    yi = pan$yi[ix],
    lci = pan$ci.lb[ix],
    uci = pan$ci.ub[ix],
    stringsAsFactors = FALSE
  )
}

forest_observed_ifrac <- function(fit) {
  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_overall = FALSE,
    include_prediction = FALSE,
    include_combined = FALSE
  )
  pan <- panels$resultsA
  labs <- panels$row_labels
  ix <- which(pan$row_kind == "study")
  data.frame(
    trial = as.character(labs$label[ix]),
    subgroup = as.character(labs$colour_label[ix]),
    ifrac = parse_ifrac_display(pan$ifrac[ix]),
    stringsAsFactors = FALSE
  )
}

expect_cams_identities <- function(cmp, tol = 1e-10) {
  for (tr in cmp$trials) {
    pi_j <- unname(cmp$pi[[tr]])
    expect_equal(
      cmp$support$y_B[, tr] - cmp$support$y_A[, tr],
      cmp$support$g[, tr],
      tolerance = tol
    )
    expect_equal(
      (1 - pi_j) * cmp$support$y_A[, tr] + pi_j * cmp$support$y_B[, tr],
      cmp$support$m[, tr],
      tolerance = tol
    )
    expect_equal(
      cmp$support$y_A[, tr],
      cmp$support$m[, tr] - pi_j * cmp$support$g[, tr],
      tolerance = tol
    )
    expect_equal(
      cmp$support$y_B[, tr],
      cmp$support$m[, tr] + (1 - pi_j) * cmp$support$g[, tr],
      tolerance = tol
    )
  }
}

expect_forest_equals_shrinkage <- function(fit) {
  sg <- shrinkage(fit, type = "modifiers", summary = TRUE, exponentiate = FALSE)
  fp <- forest_shrinkage_rows(fit)
  expect_equal(nrow(fp), nrow(sg$summary))
  key_fp <- paste(fp$trial, fp, sep = "\r")
  key_sg <- paste(sg$summary$trial, sg$summary$level, sep = "\r")
  m <- match(key_fp, key_sg)
  expect_false(any(is.na(m)))
  expect_equal(fp$yi, sg$summary$estimate[m], tolerance = 1e-8)
  expect_equal(fp$lci, sg$summary$lci[m], tolerance = 1e-6)
  expect_equal(fp$uci, sg$summary$uci[m], tolerance = 1e-6)
}

test_that("bayesmeta shrinkage equals DIRECT latent reconstruction, not FE fitted()", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  dat <- fixture_cams_native_shrinkage()
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
        silent = 2
  ))
  cmp <- shrinkage_components(fit)
  expect_equal(cmp$reference, "Younger")
  expect_equal(cmp$nonref, "Older")
  expect_cams_identities(cmp)

  sg <- shrinkage(fit, type = "modifiers", summary = TRUE, exponentiate = FALSE)
  fe <- fitted(fit, summary = TRUE)
  # Native trial latents must differ from the FE linpred at π_j
  max_fe <- 0
  for (i in seq_len(nrow(sg$summary))) {
    tr <- sg$summary$trial[[i]]
    lab <- sg$summary[[i]]
    i_dat <- which(
      as.character(fit) == tr &
        as.character(fit$data) == lab
    )[[1L]]
    max_fe <- max(max_fe, abs(sg$summary$estimate[[i]] - fe$Estimate[[i_dat]]))
  }
  expect_gt(max_fe, 0.05)

  for (tr in cmp$trials) {
    smA <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_A[[tr]])
    smB <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_B[[tr]])
    rowA <- sg$summary[
      sg$summary$trial == tr & sg$summary$level == "Younger",
      ,
      drop = FALSE
    ]
    rowB <- sg$summary[
      sg$summary$trial == tr & sg$summary$level == "Older",
      ,
      drop = FALSE
    ]
    expect_equal(rowA$estimate[[1L]], smA$median, tolerance = 1e-6)
    expect_equal(rowA$lci[[1L]], smA$lower, tolerance = 1e-5)
    expect_equal(rowA$uci[[1L]], smA$upper, tolerance = 1e-5)
    expect_equal(rowB$estimate[[1L]], smB$median, tolerance = 1e-6)
    expect_equal(rowB$lci[[1L]], smB$lower, tolerance = 1e-5)
    expect_equal(rowB$uci[[1L]], smB$upper, tolerance = 1e-5)
  }

  expect_forest_equals_shrinkage(fit)
})

test_that("metafor shrinkage equals BLUP reconstruction, not FE fitted()", {
  skip_if_not_installed("metafor")

  dat <- fixture_cams_native_shrinkage()
  fit <- suppressWarnings(cams(
    dat,
    reference = "Younger",
    method = "REML",
    silent = 2
  ))
  cmp <- shrinkage_components(fit)
  expect_cams_identities(cmp)

  sg <- shrinkage(fit, type = "modifiers", summary = TRUE, exponentiate = FALSE)
  bl_m <- metafor::blup(fit$metafor$overall)
  bl_g <- metafor::blup(fit$metafor$interaction)
  fe <- fitted(fit, summary = TRUE)
  max_fe <- 0
  for (tr in cmp$trials) {
    pi_j <- unname(cmp$pi[[tr]])
    bm <- interactmeta:::comparison_blup_row(bl_m, tr, fit$metafor$overall$slab)
    bg <- interactmeta:::comparison_blup_row(bl_g, tr, fit$metafor$interaction$slab)
    muA <- interactmeta:::interactmeta_muA_from_mg(bm$pred, bg$pred, pi_j)
    muB <- interactmeta:::interactmeta_muB_from_mg(bm$pred, bg$pred, pi_j)
    rowA <- sg$summary[
      sg$summary$trial == tr & sg$summary$level == "Younger",
      ,
      drop = FALSE
    ]
    rowB <- sg$summary[
      sg$summary$trial == tr & sg$summary$level == "Older",
      ,
      drop = FALSE
    ]
    expect_equal(rowA$estimate[[1L]], muA, tolerance = 1e-8)
    expect_equal(rowB$estimate[[1L]], muB, tolerance = 1e-8)
    iA <- which(
      as.character(fit) == tr &
        as.character(fit$data) == "Younger"
    )[[1L]]
    max_fe <- max(max_fe, abs(rowA$estimate[[1L]] - fe$Estimate[[iA]]))
  }
  expect_gt(max_fe, 0.05)
  expect_forest_equals_shrinkage(fit)
})

test_that("displayed IF is the subgroup share; reconstruction uses non-reference pi_j", {
  skip_if_not_installed("bayesmeta")

  dat <- fixture_cams_native_shrinkage()
  fit <- suppressWarnings(cams(
    dat, reference = "Younger", silent = 2
  ))
  cmp <- shrinkage_components(fit)
  disp <- forest_observed_ifrac(fit)
  for (tr in cmp$trials) {
    pi_j <- unname(cmp$pi[[tr]])
    dA <- disp$ifrac[disp$trial == tr & disp == "Younger"]
    dB <- disp$ifrac[disp$trial == tr & disp == "Older"]
    expect_equal(dB, pi_j, tolerance = 5e-4)
    expect_equal(dA, 1 - pi_j, tolerance = 5e-4)
    expect_equal(unname(cmp$pi[[tr]]), pi_j, tolerance = 1e-12)
  }
})

test_that("forest shrinkage matches shrinkage() after reordering input rows", {
  skip_if_not_installed("bayesmeta")

  dat <- fixture_cams_native_shrinkage()
  dat_rev <- fixture_cams_native_shrinkage(reverse_rows = TRUE)
  fit <- suppressWarnings(cams(
    dat, reference = "Younger", silent = 2
  ))
  fit_rev <- suppressWarnings(cams(
    dat_rev, reference = "Younger", silent = 2
  ))
  a <- shrinkage(fit, type = "modifiers", summary = TRUE, exponentiate = FALSE)$summary
  b <- shrinkage(fit_rev, type = "modifiers", summary = TRUE, exponentiate = FALSE)$summary
  key <- function(d) paste(d$trial, d, sep = "\r")
  mb <- match(key(a), key(b))
  expect_false(any(is.na(mb)))
  expect_equal(a$estimate, b$estimate[mb], tolerance = 1e-6)
  expect_forest_equals_shrinkage(fit_rev)
})

test_that("reversed reference keeps A/B reconstruction signs and forest matching", {
  skip_if_not_installed("bayesmeta")

  dat <- fixture_cams_native_shrinkage(subgroup_levels = c("Older", "Younger"))
  fit <- suppressWarnings(cams(
    dat, reference = "Older", silent = 2
  ))
  cmp <- shrinkage_components(fit)
  expect_equal(cmp$reference, "Older")
  expect_equal(cmp$nonref, "Younger")
  expect_cams_identities(cmp)
  # π_j is IF of the non-reference (Younger)
  for (tr in cmp$trials) {
    rows <- fit$data[as.character(fit) == tr, , drop = FALSE]
    if_col <- grep("^ifrac_", names(rows), value = TRUE)[[1L]]
    expect_equal(unname(cmp$pi[[tr]]), unique(rows[[if_col]])[[1L]], tolerance = 1e-10)
    i_ref <- rows == "Older"
    i_non <- rows == "Younger"
    prec <- 1 / rows$vi
    pi_non <- unname(prec[i_non] / sum(prec))
    expect_equal(unname(cmp$pi[[tr]]), pi_non, tolerance = 1e-8)
  }
  sg <- shrinkage(fit, type = "modifiers", summary = TRUE, exponentiate = FALSE)
  for (tr in cmp$trials) {
    smA <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_A[[tr]])
    rowA <- sg$summary[
      sg$summary$trial == tr & sg$summary$level == "Older",
      ,
      drop = FALSE
    ]
    expect_equal(rowA$estimate[[1L]], smA$median, tolerance = 1e-6)
  }
  expect_forest_equals_shrinkage(fit)
})

test_that("exponentiation happens only after CAMS reconstruction", {
  skip_if_not_installed("bayesmeta")

  dat <- fixture_cams_native_shrinkage()
  fit <- suppressWarnings(cams(
    dat, reference = "Younger", silent = 2
  ))
  lin <- shrinkage(fit, type = "modifiers", summary = TRUE, exponentiate = FALSE)
  expd <- shrinkage(fit, type = "modifiers", summary = TRUE, exponentiate = TRUE)
  expect_equal(expd$summary$estimate, exp(lin$summary$estimate), tolerance = 1e-10)
  expect_equal(expd$summary$lci, exp(lin$summary$lci), tolerance = 1e-10)
  expect_equal(expd$summary$uci, exp(lin$summary$uci), tolerance = 1e-10)

  cmp <- shrinkage_components(fit)
  tr <- cmp$trials[[1L]]
  pi_j <- unname(cmp$pi[[tr]])
  wrong <- exp(stats::median(cmp$support$m[, tr])) -
    pi_j * exp(stats::median(cmp$support$g[, tr]))
  ya_exp <- expd$summary$estimate[
    expd$summary$trial == tr & expd$summary$level == "Younger"
  ]
  expect_false(isTRUE(all.equal(ya_exp, wrong, tolerance = 1e-3)))

  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_overall = FALSE,
    include_prediction = FALSE,
    include_combined = FALSE
  )
  pan <- panels$resultsA
  labs <- panels$row_labels
  ix <- which(
    pan$row_kind == "shrinkage" &
      labs$label == tr &
      labs$colour_label == "Younger" &
      !is.na(pan$shape)
  )
  expect_equal(pan$yi[[ix[[1L]]]], lin$summary$estimate[
    lin$summary$trial == tr & lin$summary$level == "Younger"
  ][[1L]], tolerance = 1e-8)
})
