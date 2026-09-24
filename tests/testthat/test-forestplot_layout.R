test_that("forestplot_figure_inches scales with rows and footnotes", {
  fp <- structure(
    list(
      engine = "metafor",
      specs = list(
        panels = list(data.frame(yi = seq_len(10L))),
        row_labels = data.frame(
          heading = rep(c("A", "B"), each = 5L),
          label = as.character(seq_len(10L)),
          stringsAsFactors = FALSE
        ),
        annotations = list(
          lines = character(),
          ifrac = character(),
          heterogeneity = character()
        )
      ),
      figure = structure(list(engine = "metafor"), class = "interactmeta_forestplot_metafor")
    ),
    class = "interactmeta_forestplot",
    panels = list(data.frame(yi = seq_len(10L)))
  )
  dims <- interactmeta:::forestplot_figure_inches(fp)
  # Full forests: ~0.28 in per layout row + ~2 in chrome.
  expect_match(dims$formula, "0\\.28\\*")
  expect_false(grepl("aspect floor", dims$formula, fixed = TRUE))
  expect_gte(dims$width, 12)
  h <- forestplot_height(fp, width = 8)
  expect_equal(h$width, 8)
  expect_gte(h$height, 3)
  expect_lt(h$height, 0.9 * 12.5)
  expect_true(is.finite(h$n_rows) && h$n_rows >= 10L)

  # A caller-chosen width must not inherit height inflated from the default
  # 12.5 in recommendation (that made vignette forests too tall).
  h10 <- forestplot_height(fp, width = 10)
  expect_equal(h10$width, 10)
  expect_equal(h10$height, h$height)

  # Annotation footnotes are counted for sizing.
  fp$specs$annotations <- list(
    lines = character(),
    ifrac = c(Older = "Information fraction (\u03c0): 0.40 (0.30, 0.50)"),
    heterogeneity = c(
      overall = "Heterogeneity (\u03c4): 0.10 (0.01, 0.30)",
      interaction = "Interaction heterogeneity (\u03c4_\u03b3): 0.05 (0.01, 0.20)"
    )
  )
  dims_ann <- interactmeta:::forestplot_figure_inches(fp)
  expect_gte(dims_ann$n_footnotes, 1L)
  expect_gt(dims_ann$n_footnotes, dims$n_footnotes)
})

test_that("forestplot_height scales for a large multi-subgroup forest", {
  # ~4 blocks × 7 studies + pooled/prediction ≈ age-like default forest.
  n <- 36L
  headings <- rep(c("Younger", "Older", "Interaction", "Overall"), each = 9L)
  fp <- structure(
    list(
      engine = "metafor",
      specs = list(
        panels = list(data.frame(
          yi = seq_len(n),
          ci.lb = seq_len(n) - 0.2,
          ci.ub = seq_len(n) + 0.2,
          row_kind = rep(c(rep("study", 7L), "pooled", "prediction"), 4L)
        )),
        row_labels = data.frame(
          heading = headings,
          label = paste0("row", seq_len(n)),
          stringsAsFactors = FALSE
        ),
        estimate_col = "yi"
      ),
      figure = structure(list(engine = "metafor"), class = "interactmeta_forestplot_metafor")
    ),
    class = "interactmeta_forestplot",
    panels = list(data.frame(
      yi = seq_len(n),
      ci.lb = seq_len(n) - 0.2,
      ci.ub = seq_len(n) + 0.2,
      row_kind = rep(c(rep("study", 7L), "pooled", "prediction"), 4L)
    ))
  )
  dims <- forestplot_height(fp, width = 10)
  expect_equal(dims$width, 10)
  expect_gte(dims$n_rows, n)
  expect_true(!is.null(dims$formula))
  expect_gte(dims$height, 0.28 * (dims$n_rows + 1L) + 2)
  expect_lte(dims$height, 24)
  # Fixed fig.height = 6 would be far too short for this forest.
  expect_gt(dims$height, 8)
})
test_that("forestplot_metafor_cex() shrinks on a short device", {
  tall <- interactmeta:::forestplot_metafor_cex(9, n_layout_rows = 20, device_height_in = 10)
  short <- interactmeta:::forestplot_metafor_cex(9, n_layout_rows = 20, device_height_in = 3)
  expect_lt(short, tall)
  expect_lte(short, 0.75)
})

test_that("forestplot_display_metrics() scales with device and rows", {
  m_tall <- interactmeta:::forestplot_display_metrics(
    n_rows = 20,
    device_width_in = 12,
    device_height_in = 10,
    has_favours = TRUE
  )
  m_short <- interactmeta:::forestplot_display_metrics(
    n_rows = 20,
    device_width_in = 12,
    device_height_in = 4,
    has_favours = TRUE
  )
  expect_lt(m_short$base_size, m_tall$base_size)
  expect_true(is.finite(m_tall$panel.height))
  expect_true(is.finite(m_short$panel.height))
})

test_that("forestplot_n_layout_rows counts headers", {
  fp <- structure(
    list(
      engine = "metafor",
      specs = list(
        panels = list(data.frame(yi = 1:4)),
        row_labels = data.frame(
          heading = c("A", "A", "B", "B"),
          label = letters[1:4],
          stringsAsFactors = FALSE
        )
      )
    ),
    class = "interactmeta_forestplot",
    panels = list(data.frame(yi = 1:4))
  )
  # 4 content + 2 unique headings + 1 table header
  expect_equal(interactmeta:::forestplot_n_layout_rows(fp), 7L)
})
