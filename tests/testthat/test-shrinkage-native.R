# Engine-native CAMS shrinkage (bayesmeta DIRECT / metafor BLUP).
# Public engines have no r_trial columns.

cams_native_shrinkage_data <- function(reverse_rows = FALSE,
                                       subgroup_levels = c("Younger", "Older")) {
  dat <- data.frame(
    trial = rep(c("ADVANCE", "AFFIRM", "DEFINE", "TEMSO"), each = 2),
    subgroup = factor(
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
    include_pooled = FALSE,
    include_prediction = FALSE,
    include_overall = FALSE,
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
    include_pooled = FALSE,
    include_prediction = FALSE,
    include_overall = FALSE
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
      cmp$draws$y_B[, tr] - cmp$draws$y_A[, tr],
      cmp$draws$g[, tr],
      tolerance = tol
    )
    expect_equal(
      (1 - pi_j) * cmp$draws$y_A[, tr] + pi_j * cmp$draws$y_B[, tr],
      cmp$draws$m[, tr],
      tolerance = tol
    )
    expect_equal(
      cmp$draws$y_A[, tr],
      cmp$draws$m[, tr] - pi_j * cmp$draws$g[, tr],
      tolerance = tol
    )
    expect_equal(
      cmp$draws$y_B[, tr],
      cmp$draws$m[, tr] + (1 - pi_j) * cmp$draws$g[, tr],
      tolerance = tol
    )
  }
}

expect_forest_equals_shrinkage <- function(fit) {
  sg <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = FALSE)
  fp <- forest_shrinkage_rows(fit)
  expect_equal(nrow(fp), nrow(sg$summary))
  key_fp <- paste(fp$trial, fp$subgroup, sep = "\r")
  key_sg <- paste(sg$summary$trial, sg$summary$subgroup, sep = "\r")
  m <- match(key_fp, key_sg)
  expect_false(any(is.na(m)))
  expect_equal(fp$yi, sg$summary$estimate[m], tolerance = 1e-8)
  expect_equal(fp$lci, sg$summary$lci[m], tolerance = 1e-6)
  expect_equal(fp$uci, sg$summary$uci[m], tolerance = 1e-6)
}

test_that("bayesmeta shrinkage equals DIRECT latent reconstruction, not FE fitted()", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("metafor")

  dat <- cams_native_shrinkage_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "bayesmeta",
    silent = 2
  ))
  cmp <- shrinkage_components(fit)
  expect_equal(cmp$reference, "Younger")
  expect_equal(cmp$nonref, "Older")
  expect_cams_identities(cmp)

  sg <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = FALSE)
  fe <- fitted(fit, summary = TRUE)
  # Native trial latents must differ from the FE linpred at π_j
  max_fe <- 0
  for (i in seq_len(nrow(sg$summary))) {
    tr <- sg$summary$trial[[i]]
    lab <- sg$summary$subgroup[[i]]
    i_dat <- which(
      as.character(fit$data$trial) == tr &
        as.character(fit$data$subgroup) == lab
    )[[1L]]
    max_fe <- max(max_fe, abs(sg$summary$estimate[[i]] - fe$Estimate[[i_dat]]))
  }
  expect_gt(max_fe, 0.05)

  for (tr in cmp$trials) {
    smA <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_A[[tr]])
    smB <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_B[[tr]])
    rowA <- sg$summary[
      sg$summary$trial == tr & sg$summary$subgroup == "Younger",
      ,
      drop = FALSE
    ]
    rowB <- sg$summary[
      sg$summary$trial == tr & sg$summary$subgroup == "Older",
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

  dat <- cams_native_shrinkage_data()
  fit <- suppressWarnings(cam(
    dat,
    reference = "Younger",
    engine = "metafor",
    silent = 2
  ))
  cmp <- shrinkage_components(fit)
  expect_cams_identities(cmp)

  sg <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = FALSE)
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
      sg$summary$trial == tr & sg$summary$subgroup == "Younger",
      ,
      drop = FALSE
    ]
    rowB <- sg$summary[
      sg$summary$trial == tr & sg$summary$subgroup == "Older",
      ,
      drop = FALSE
    ]
    expect_equal(rowA$estimate[[1L]], muA, tolerance = 1e-8)
    expect_equal(rowB$estimate[[1L]], muB, tolerance = 1e-8)
    iA <- which(
      as.character(fit$data$trial) == tr &
        as.character(fit$data$subgroup) == "Younger"
    )[[1L]]
    max_fe <- max(max_fe, abs(rowA$estimate[[1L]] - fe$Estimate[[iA]]))
  }
  expect_gt(max_fe, 0.05)
  expect_forest_equals_shrinkage(fit)
})

test_that("displayed IF is the subgroup share; reconstruction uses non-reference pi_j", {
  skip_if_not_installed("bayesmeta")

  dat <- cams_native_shrinkage_data()
  fit <- suppressWarnings(cam(
    dat, reference = "Younger", engine = "bayesmeta", silent = 2
  ))
  cmp <- shrinkage_components(fit)
  disp <- forest_observed_ifrac(fit)
  for (tr in cmp$trials) {
    pi_j <- unname(cmp$pi[[tr]])
    dA <- disp$ifrac[disp$trial == tr & disp$subgroup == "Younger"]
    dB <- disp$ifrac[disp$trial == tr & disp$subgroup == "Older"]
    expect_equal(dB, pi_j, tolerance = 5e-4)
    expect_equal(dA, 1 - pi_j, tolerance = 5e-4)
    expect_equal(unname(cmp$pi[[tr]]), pi_j, tolerance = 1e-12)
  }
})

test_that("forest shrinkage matches shrinkage() after reordering input rows", {
  skip_if_not_installed("bayesmeta")

  dat <- cams_native_shrinkage_data()
  dat_rev <- cams_native_shrinkage_data(reverse_rows = TRUE)
  fit <- suppressWarnings(cam(
    dat, reference = "Younger", engine = "bayesmeta", silent = 2
  ))
  fit_rev <- suppressWarnings(cam(
    dat_rev, reference = "Younger", engine = "bayesmeta", silent = 2
  ))
  a <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = FALSE)$summary
  b <- shrinkage(fit_rev, type = "subgroups", summary = TRUE, exponentiate = FALSE)$summary
  key <- function(d) paste(d$trial, d$subgroup, sep = "\r")
  mb <- match(key(a), key(b))
  expect_false(any(is.na(mb)))
  expect_equal(a$estimate, b$estimate[mb], tolerance = 1e-6)
  expect_forest_equals_shrinkage(fit_rev)
})

test_that("reversed reference keeps A/B reconstruction signs and forest matching", {
  skip_if_not_installed("bayesmeta")

  dat <- cams_native_shrinkage_data(subgroup_levels = c("Older", "Younger"))
  fit <- suppressWarnings(cam(
    dat, reference = "Older", engine = "bayesmeta", silent = 2
  ))
  cmp <- shrinkage_components(fit)
  expect_equal(cmp$reference, "Older")
  expect_equal(cmp$nonref, "Younger")
  expect_cams_identities(cmp)
  # π_j is IF of the non-reference (Younger)
  for (tr in cmp$trials) {
    rows <- fit$data[as.character(fit$data$trial) == tr, , drop = FALSE]
    if_col <- grep("^ifrac_", names(rows), value = TRUE)[[1L]]
    expect_equal(unname(cmp$pi[[tr]]), unique(rows[[if_col]])[[1L]], tolerance = 1e-10)
    i_ref <- rows$subgroup == "Older"
    i_non <- rows$subgroup == "Younger"
    prec <- 1 / rows$vi
    pi_non <- unname(prec[i_non] / sum(prec))
    expect_equal(unname(cmp$pi[[tr]]), pi_non, tolerance = 1e-8)
  }
  sg <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = FALSE)
  for (tr in cmp$trials) {
    smA <- interactmeta:::interactmeta_posterior_1d_summary(fit$latent$mu_A[[tr]])
    rowA <- sg$summary[
      sg$summary$trial == tr & sg$summary$subgroup == "Older",
      ,
      drop = FALSE
    ]
    expect_equal(rowA$estimate[[1L]], smA$median, tolerance = 1e-6)
  }
  expect_forest_equals_shrinkage(fit)
})

test_that("exponentiation happens only after CAMS reconstruction", {
  skip_if_not_installed("bayesmeta")

  dat <- cams_native_shrinkage_data()
  fit <- suppressWarnings(cam(
    dat, reference = "Younger", engine = "bayesmeta", silent = 2
  ))
  lin <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = FALSE)
  expd <- shrinkage(fit, type = "subgroups", summary = TRUE, exponentiate = TRUE)
  expect_equal(expd$summary$estimate, exp(lin$summary$estimate), tolerance = 1e-10)
  expect_equal(expd$summary$lci, exp(lin$summary$lci), tolerance = 1e-10)
  expect_equal(expd$summary$uci, exp(lin$summary$uci), tolerance = 1e-10)

  cmp <- shrinkage_components(fit)
  tr <- cmp$trials[[1L]]
  pi_j <- unname(cmp$pi[[tr]])
  wrong <- exp(stats::median(cmp$draws$m[, tr])) -
    pi_j * exp(stats::median(cmp$draws$g[, tr]))
  ya_exp <- expd$summary$estimate[
    expd$summary$trial == tr & expd$summary$subgroup == "Younger"
  ]
  expect_false(isTRUE(all.equal(ya_exp, wrong, tolerance = 1e-3)))

  panels <- interactmeta:::comparison_forest_build_fit_panels(
    fit,
    prob = 0.95,
    contrast_sep = "-",
    max_ifrac_draws = NULL,
    model_names = NULL,
    include_shrinkage = TRUE,
    include_pooled = FALSE,
    include_prediction = FALSE,
    include_overall = FALSE
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
    lin$summary$trial == tr & lin$summary$subgroup == "Younger"
  ][[1L]], tolerance = 1e-8)
})
