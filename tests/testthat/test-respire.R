respire_published_age <- function() {
  data.frame(
    trial = c("RESPIRE 1", "RESPIRE 1", "RESPIRE 2", "RESPIRE 2"),
    subgroup = c("<65", ">=65", "<65", ">=65"),
    n = c(107L, 168L, 200L, 150L),
    hr = c(0.61, 0.49, 0.69, 1.19),
    ci_lower = c(0.32, 0.29, 0.45, 0.70),
    ci_upper = c(1.16, 0.84, 1.07, 2.03),
    ci_level = c(0.975, 0.975, 0.951, 0.951),
    stringsAsFactors = FALSE
  )
}

respire_published_sex <- function() {
  data.frame(
    trial = c("RESPIRE 1", "RESPIRE 1", "RESPIRE 2", "RESPIRE 2"),
    subgroup = c("Female", "Male", "Female", "Male"),
    n = c(184L, 91L, 210L, 140L),
    hr = c(0.49, 0.63, 0.64, 1.49),
    ci_lower = c(0.30, 0.31, 0.40, 0.84),
    ci_upper = c(0.81, 1.29, 1.00, 2.63),
    ci_level = c(0.975, 0.975, 0.951, 0.951),
    stringsAsFactors = FALSE
  )
}

respire_sort <- function(d) {
  d[order(as.character(d$trial), as.character(d$subgroup)), , drop = FALSE]
}

test_that("respire_age retains published values and derived yi/sei", {
  src <- respire_published_age()
  expect_equal(nrow(respire_age), 4L)
  expect_equal(as.character(respire_age$trial), src$trial)
  expect_equal(as.character(respire_age$subgroup), src$subgroup)
  expect_equal(respire_age$n, src$n)
  expect_equal(respire_age$hr, src$hr)
  expect_equal(respire_age$ci_lower, src$ci_lower)
  expect_equal(respire_age$ci_upper, src$ci_upper)
  expect_equal(respire_age$ci_level, src$ci_level)
  expect_equal(levels(respire_age$subgroup), c("<65", ">=65"))

  trans <- log_hr_from_ci(
    respire_age$hr,
    respire_age$ci_lower,
    respire_age$ci_upper,
    respire_age$ci_level
  )
  expect_equal(respire_age$yi, log(respire_age$hr), tolerance = 1e-12)
  expect_equal(respire_age$yi, trans$yi, tolerance = 1e-12)
  expect_equal(respire_age$sei, trans$sei, tolerance = 1e-12)
  expect_equal(respire_age$vi, trans$vi, tolerance = 1e-12)

  recon <- hr_ci_from_log(respire_age$yi, respire_age$sei, respire_age$ci_level)
  expect_equal(recon$ci_lower, respire_age$ci_lower, tolerance = 0.02)
  expect_equal(recon$ci_upper, respire_age$ci_upper, tolerance = 0.02)

  z_wrong <- 1.96
  sei_wrong <- (log(respire_age$ci_upper) - log(respire_age$ci_lower)) / (2 * z_wrong)
  expect_false(isTRUE(all.equal(sei_wrong, respire_age$sei, tolerance = 1e-8)))
})

test_that("respire_sex retains published values and derived yi/sei", {
  src <- respire_published_sex()
  expect_equal(nrow(respire_sex), 4L)
  expect_equal(as.character(respire_sex$trial), src$trial)
  expect_equal(as.character(respire_sex$subgroup), src$subgroup)
  expect_equal(respire_sex$n, src$n)
  expect_equal(respire_sex$hr, src$hr)
  expect_equal(respire_sex$ci_lower, src$ci_lower)
  expect_equal(respire_sex$ci_upper, src$ci_upper)
  expect_equal(respire_sex$ci_level, src$ci_level)
  expect_equal(levels(respire_sex$subgroup), c("Male", "Female"))

  trans <- log_hr_from_ci(
    respire_sex$hr,
    respire_sex$ci_lower,
    respire_sex$ci_upper,
    respire_sex$ci_level
  )
  expect_equal(respire_sex$yi, log(respire_sex$hr), tolerance = 1e-12)
  expect_equal(respire_sex$sei, trans$sei, tolerance = 1e-12)
  recon <- hr_ci_from_log(respire_sex$yi, respire_sex$sei, respire_sex$ci_level)
  expect_equal(recon$ci_lower, respire_sex$ci_lower, tolerance = 0.02)
  expect_equal(recon$ci_upper, respire_sex$ci_upper, tolerance = 0.02)
})

test_that("RESPIRE 1 uses 97.5% intervals and RESPIRE 2 uses 95.1%", {
  expect_equal(unique(respire_age$ci_level[respire_age$trial == "RESPIRE 1"]), 0.975)
  expect_equal(unique(respire_age$ci_level[respire_age$trial == "RESPIRE 2"]), 0.951)
  expect_equal(unique(respire_sex$ci_level[respire_sex$trial == "RESPIRE 1"]), 0.975)
  expect_equal(unique(respire_sex$ci_level[respire_sex$trial == "RESPIRE 2"]), 0.951)
})

test_that("manual log-HR/SE reconstruction matches packaged respire_age", {
  x <- subset(respire_age, trial == "RESPIRE 1" & as.character(subgroup) == "<65")
  yi <- log(x$hr)
  z <- qnorm(1 - (1 - x$ci_level) / 2)
  sei <- (log(x$ci_upper) - log(x$ci_lower)) / (2 * z)
  expect_equal(as.numeric(x$yi), as.numeric(yi), tolerance = 1e-12)
  expect_equal(as.numeric(x$sei), as.numeric(sei), tolerance = 1e-12)
})

fit_respire_three <- function(data, reference) {
  skip_if_not_installed("metafor")
  skip_if_not_installed("bayesmeta")
  conv <- conventional(data, reference = reference)
  sw <- swada(data, reference = reference)
  cam <- suppressWarnings(cam(
    data,
    reference = reference,
    engine = "bayesmeta",
    silent = 2
  ))
  list(conventional = conv, swada = sw, cam = cam)
}

test_that("all three age approaches use the same published yi/vi input", {
  fits <- fit_respire_three(respire_age, "<65")
  src <- respire_sort(respire_age[, c("trial", "subgroup", "yi", "vi")])
  conv_in <- respire_sort(fits$conventional$data[, c("trial", "subgroup", "yi", "vi")])
  sw_in <- respire_sort(fits$swada$data[, c("trial", "subgroup", "yi", "vi")])
  cams_in <- respire_sort(fits$cam$data[, c("trial", "subgroup", "yi", "vi")])
  expect_equal(conv_in$yi, src$yi, tolerance = 1e-12)
  expect_equal(sw_in$yi, src$yi, tolerance = 1e-12)
  expect_equal(cams_in$yi, src$yi, tolerance = 1e-12)
  expect_equal(conv_in$vi, src$vi, tolerance = 1e-12)
  expect_equal(sw_in$vi, src$vi, tolerance = 1e-12)
  expect_equal(cams_in$vi, src$vi, tolerance = 1e-12)
  expect_equal(fits$conventional$reference, "<65")
  expect_equal(fits$swada$reference, "<65")
  expect_equal(fits$cam$model_spec$reference, "<65")
})

test_that("all three sex approaches use the same published yi/vi input", {
  fits <- fit_respire_three(respire_sex, "Male")
  src <- respire_sort(respire_sex[, c("trial", "subgroup", "yi", "vi")])
  conv_in <- respire_sort(fits$conventional$data[, c("trial", "subgroup", "yi", "vi")])
  sw_in <- respire_sort(fits$swada$data[, c("trial", "subgroup", "yi", "vi")])
  cams_in <- respire_sort(fits$cam$data[, c("trial", "subgroup", "yi", "vi")])
  expect_equal(conv_in$yi, src$yi, tolerance = 1e-12)
  expect_equal(sw_in$yi, src$yi, tolerance = 1e-12)
  expect_equal(cams_in$yi, src$yi, tolerance = 1e-12)
  expect_equal(fits$conventional$reference, "Male")
  expect_equal(fits$swada$reference, "Male")
  expect_equal(fits$cam$model_spec$reference, "Male")
})

test_that("CAMS age transformation matches g = y_>=65 - y_<65 and published SEs", {
  fits <- fit_respire_three(respire_age, "<65")
  q <- cam_trial_quantities(fits$cam)
  pairs <- fits$conventional$pairs
  expect_equal(attr(q, "reference"), "<65")
  expect_equal(attr(q, "nonref"), ">=65")
  expect_equal(q$g, pairs$g, tolerance = 1e-12)
  expect_equal(q$pi, pairs$pi, tolerance = 1e-12)
  expect_equal(q$m, pairs$m, tolerance = 1e-12)
  r1 <- subset(respire_age, trial == "RESPIRE 1")
  g_manual <- as.numeric(r1$yi[r1$subgroup == ">=65"] - r1$yi[r1$subgroup == "<65"])
  expect_equal(q$g[q$trial == "RESPIRE 1"], g_manual, tolerance = 1e-12)
})

test_that("CAMS sex transformation matches g = y_Female - y_Male", {
  fits <- fit_respire_three(respire_sex, "Male")
  q <- cam_trial_quantities(fits$cam)
  expect_equal(attr(q, "reference"), "Male")
  expect_equal(attr(q, "nonref"), "Female")
  r1 <- subset(respire_sex, trial == "RESPIRE 1")
  g_manual <- as.numeric(r1$yi[r1$subgroup == "Female"] - r1$yi[r1$subgroup == "Male"])
  expect_equal(q$g[q$trial == "RESPIRE 1"], g_manual, tolerance = 1e-12)
})

test_that("SWADA uses one weight vector; DA weights differ by subgroup series", {
  fits <- fit_respire_three(respire_age, "<65")
  w <- compare_weights(fits$conventional, fits$swada)
  expect_equal(nrow(w), 2L)
  expect_equal(sum(w$swada), 1, tolerance = 1e-12)
  expect_equal(w$swada, w$ad, tolerance = 1e-10)
  expect_true(isTRUE(fits$swada$collapsibility$equal))
})

test_that("comparison tables use DA, AD, SWADA, and CAMS labels", {
  fits <- fit_respire_three(respire_age, "<65")
  tab <- suppressWarnings(compare_subgroups(fits$conventional, fits$swada, fits$cam))
  expect_equal(
    sort(unique(as.character(tab$approach))),
    c("CAMS", "DA", "SWADA")
  )
  expect_equal(sort(unique(as.character(tab$subgroup))), c("<65", ">=65"))
  expect_equal(nrow(tab), 6L)
  expect_true(all(is.finite(tab$estimate)))
  ix <- suppressWarnings(compare_interactions(fits$conventional, fits$swada, fits$cam))
  expect_equal(nrow(ix), 4L)
  expect_equal(
    as.character(ix$approach),
    c("DA", "AD", "SWADA", "CAMS")
  )
  expect_true(all(is.finite(ix$estimate)))
})

test_that("package datasets include respire_age and respire_sex", {
  items <- utils::data(package = "interact")$results[, "Item"]
  expect_true(all(c("respire_age", "respire_sex") %in% items))
})
