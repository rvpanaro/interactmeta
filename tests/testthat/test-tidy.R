test_that("tidy.interactmeta_subgroup_effects() returns summary columns", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("generics")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger")
  me <- subgroups(fit)
  td <- generics::tidy(me)

  expect_s3_class(td, "tbl_df")
  expect_true(all(c("estimate", "std.error", "conf.low", "conf.high") %in% names(td)))
  expect_true(all(c("subgroup", "type") %in% names(td)))
  expect_false("term" %in% names(td))
  expect_false("lci" %in% names(td))
  expect_false("sd" %in% names(td))
  expect_equal(nrow(td), nrow(me$summary))
})

test_that("tidy.interactmeta_interaction_effects() deduplicates shrinkage summaries", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("generics")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(dat, engine = "bayesmeta", reference = "Younger")
  me <- interactions(fit, shrinkage = TRUE)
  td <- generics::tidy(me)

  expect_equal(nrow(td), 1L)
  expect_equal(unique(td$interaction), "Older-Younger")
})

test_that("print methods use biostatistical tables (not tibbles)", {
  skip_if_not_installed("bayesmeta")
  skip_if_not_installed("generics")

  dat <- data.frame(
    trial = rep(c("A", "B"), each = 2),
    subgroup = factor(rep(c("Younger", "Older"), 2), levels = c("Younger", "Older")),
    yi = c(0.1, 0.2, 0.15, 0.25),
    vi = 0.01,
    reference = "Younger",
    ifrac_Older = c(0.3, 0.3, 0.4, 0.4),
    ifrac2_Older = c(-0.3, 0.7, -0.4, 0.6)
  )
  fit <- cam(
    dat,
    engine = "bayesmeta",
    reference = "Younger")

  out <- capture.output(print(subgroups(fit)))
  expect_true(any(grepl("subgroup", out, ignore.case = TRUE)))
  expect_true(any(grepl("95% CrI", out, fixed = TRUE)))
  expect_true(any(grepl("95% Pred", out, fixed = TRUE)))
  expect_false(any(grepl("# A tibble", out, fixed = TRUE)))
  expect_false(any(grepl("^\\s*Prediction\\s+", out)))
  expect_false(any(grepl("^\\s*term\\b", out)))

  ov <- capture.output(print(overall(fit)))
  expect_true(any(grepl("95% Pred", ov, fixed = TRUE)))
  expect_false(any(grepl("^\\s*Prediction\\s+", ov)))
  expect_false(any(grepl("^\\s*Pooled\\s+", ov)))
})
