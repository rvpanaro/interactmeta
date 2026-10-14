test_that("interactmeta_palette() returns named colours (Okabe–Ito by default)", {
  pal <- interactmeta_palette()
  expect_named(
    pal,
    c("teal", "navy", "teal_light", "mint", "navy_mid", "slate", "white")
  )
  expect_match(pal[["teal"]], "^#[0-9A-Fa-f]{6}$")
  expect_equal(pal[["white"]], "#FFFFFF")
  expect_equal(pal[["slate"]], "#6B7C8A")
  expect_equal(pal, interactmeta_palette(palette = "okabeito"))
})

test_that("subgroup_colour_map() assigns reference and non-reference colours", {
  cols <- interact:::subgroup_colour_map(c("Lower", "Higher", "Lower"), reference = "Lower")
  by_level <- interact:::subgroup_colour_map(c("Lower", "Higher"), reference = "Lower")
  expect_equal(cols$colour[1], by_level$colour[1])
  expect_equal(cols$colour[2], by_level$colour[2])
  expect_equal(cols$fill[1], by_level$fill[1])
  expect_equal(cols$colour[1], cols$colour[3])
  # Default: red then blue (reference / non-reference)
  expect_equal(interact:::interactmeta_normalize_hex(by_level$colour[1]), "#E41A1C")
  expect_equal(interact:::interactmeta_normalize_hex(by_level$colour[2]), "#377EB8")
})

test_that("interaction_row_colours() uses non-reference subgroup colours only", {
  skip_if_not_installed("paletteer")
  by_level <- interact:::subgroup_colour_map(c("Lower", "Higher"), reference = "Lower")
  higher_line <- by_level$colour[2]
  cols <- interact:::interaction_row_colours(
    labels = c("Lower", "Higher", "Higher-Lower"),
    yi = c(NA, 0.1, 0.1),
    reference = "Lower",
    contrast_sep = "-"
  )
  expect_true(is.na(cols$colour[1]))
  expect_equal(cols$colour[2], higher_line)
  expect_equal(cols$colour[3], higher_line)
})

test_that("subgroup_shrinkage_colour_map() uses lighter grey for reference", {
  skip_if_not_installed("paletteer")
  cols <- interact:::subgroup_shrinkage_colour_map(
    c("Younger", "Older"),
    reference = "Younger"
  )
  expect_false(identical(cols$fill[1], cols$fill[2]))
  expect_match(cols$fill[1], "^#[0-9A-Fa-f]{6}$")
  expect_match(cols$colour[1], "^#[0-9A-Fa-f]{6}$")
  expect_true(as.numeric(strtoi(substring(cols$colour[1], 2, 3), 16L)) <
                as.numeric(strtoi(substring(cols$fill[1], 2, 3), 16L)))
})

test_that("interaction_row_colours() skips empty labels without error", {
  skip_if_not_installed("paletteer")
  by_level <- interact:::subgroup_colour_map(c("Younger", "Older"), reference = "Younger")
  other_line <- by_level$colour[2]
  cols <- interact:::interaction_row_colours(
    labels = c("Younger", "", "Older"),
    yi = c(0.1, 0.2, 0.15),
    reference = "Younger"
  )
  expect_equal(cols$colour[2], NA_character_)
  expect_equal(cols$colour[3], other_line)
})

test_that("subgroupmeta_palette() remains an alias of interactmeta_palette()", {
  expect_identical(subgroupmeta_palette(), interactmeta_palette())
})

test_that("plot chrome does not reuse subgroup identity colours", {
  const <- interact:::interactmeta_plot_constants()
  pal <- interactmeta_palette()
  sg <- interact:::interactmeta_okabe_ito_colours(2L)
  het <- interact:::plot_heterogeneity_colours()
  ap <- interact:::interactmeta_approach_colours()
  expect_equal(const$forest_summary, unname(pal[["navy"]]))
  expect_equal(const$forest_study, "#000000")
  expect_false(unname(het[["overall"]]) %in% sg)
  expect_false(unname(het[["interaction"]]) %in% sg)
  expect_equal(unname(ap[["Published"]]), const$text_muted)
  expect_equal(unname(ap[["DA"]]), unname(pal[["navy"]]))
  expect_equal(unname(ap[["SWADA"]]), unname(pal[["teal"]]))
  expect_equal(unname(ap[["CAMS"]]), unname(pal[["navy_mid"]]))
})
