test_that("CAMS identities hold by independent calculation (equal variances)", {
  yA <- 0.10
  yB <- 0.30
  vA <- 0.04
  vB <- 0.04
  pi <- (1 / vB) / (1 / vA + 1 / vB)
  m <- (1 - pi) * yA + pi * yB
  g <- yB - yA
  vm <- 1 / (1 / vA + 1 / vB)
  vg <- vA + vB

  expect_equal(pi, 0.5, tolerance = 1e-12)
  expect_equal(m, 0.20, tolerance = 1e-12)
  expect_equal(g, 0.20, tolerance = 1e-12)
  expect_equal(vm, 0.02, tolerance = 1e-12)
  expect_equal(vg, 0.08, tolerance = 1e-12)
  expect_equal(-(1 - pi) * vA + pi * vB, 0, tolerance = 1e-12)

  muA <- m - pi * g
  muB <- m + (1 - pi) * g
  expect_equal(muA, yA, tolerance = 1e-12)
  expect_equal(muB, yB, tolerance = 1e-12)

  dat <- data.frame(
    trial = "T",
    subgroup = c("A", "B"),
    yi = c(yA, yB),
    vi = c(vA, vB),
    ifrac_B = pi
  )
  blocks <- interact:::comparison_trial_effect_blocks(
    data = dat,
    reference = "A",
    nonref = "B",
    ifrac_var = "ifrac_B"
  )
  expect_equal(blocks$overall$yi, m, tolerance = 1e-12)
  expect_equal(blocks$overall$vi, vm, tolerance = 1e-12)
  expect_equal(blocks$interaction$yi, g, tolerance = 1e-12)
  expect_equal(blocks$interaction$vi, vg, tolerance = 1e-12)
  expect_equal(blocks$overall$ifrac, pi, tolerance = 1e-12)
})

test_that("CAMS identities hold with unequal variances and unbalanced IF", {
  yA <- 0.05
  yB <- 0.40
  vA <- 0.09
  vB <- 0.01
  pi <- (1 / vB) / (1 / vA + 1 / vB)
  m <- (1 - pi) * yA + pi * yB
  g <- yB - yA
  vm <- 1 / (1 / vA + 1 / vB)
  vg <- vA + vB

  expect_equal(pi, 0.9, tolerance = 1e-12)
  expect_equal(-(1 - pi) * vA + pi * vB, 0, tolerance = 1e-12)
  expect_equal(m - pi * g, yA, tolerance = 1e-12)
  expect_equal(m + (1 - pi) * g, yB, tolerance = 1e-12)
  expect_equal(g, 0.35, tolerance = 1e-12)

  dat <- data.frame(
    trial = "T",
    subgroup = c("A", "B"),
    yi = c(yA, yB),
    vi = c(vA, vB),
    ifrac_B = pi
  )
  blocks <- interact:::comparison_trial_effect_blocks(
    data = dat,
    reference = "A",
    nonref = "B",
    ifrac_var = "ifrac_B"
  )
  expect_equal(unname(blocks$overall$yi), m, tolerance = 1e-12)
  expect_equal(unname(blocks$overall$vi), vm, tolerance = 1e-12)
  expect_equal(unname(blocks$interaction$yi), g, tolerance = 1e-12)
  expect_equal(unname(blocks$interaction$vi), vg, tolerance = 1e-12)
})

test_that("zero interaction and balanced IF recover a common mean", {
  yA <- 0.22
  yB <- 0.22
  vA <- 0.05
  vB <- 0.05
  pi <- 0.5
  m <- 0.22
  g <- 0
  expect_equal((1 - pi) * yA + pi * yB, m)
  expect_equal(yB - yA, g)
  expect_equal(m - pi * g, yA)
  expect_equal(m + (1 - pi) * g, yB)
})

test_that("nonzero interaction with equal variances is a location shift", {
  yA <- -0.10
  yB <- 0.30
  vA <- 0.02
  vB <- 0.02
  pi <- 0.5
  g <- 0.40
  m <- 0.10
  expect_equal(yB - yA, g)
  expect_equal((1 - pi) * yA + pi * yB, m)
  expect_equal(m - pi * g, yA)
  expect_equal(m + (1 - pi) * g, yB)
  expect_equal(-(1 - pi) * vA + pi * vB, 0)
})

test_that("slope identity theta = beta + gamma is algebraic", {
  alpha <- 0.1
  beta <- -0.2
  gamma <- 0.5
  theta <- beta + gamma
  expect_equal(theta - gamma, beta)
  expect_equal(alpha + theta * 0.3, alpha + beta * 0.3 + gamma * 0.3)
})
