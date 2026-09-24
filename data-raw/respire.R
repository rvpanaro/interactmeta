# Published RESPIRE 1 / RESPIRE 2 subgroup HRs for time to first
# exacerbation, ciprofloxacin DPI 14-day on/off vs pooled placebo.
# Keep the original published HR / CI / CI-level columns; derive yi/sei
# with the same formula as log_hr_from_ci() (do not use z = 1.96).

respire_derive <- function(hr, ci_lower, ci_upper, ci_level) {
  alpha <- 1 - ci_level
  z <- stats::qnorm(1 - alpha / 2)
  yi <- log(hr)
  sei <- (log(ci_upper) - log(ci_lower)) / (2 * z)
  list(yi = yi, sei = sei, vi = sei^2)
}

respire_age_src <- data.frame(
  trial = c("RESPIRE 1", "RESPIRE 1", "RESPIRE 2", "RESPIRE 2"),
  subgroup = factor(c("<65", ">=65", "<65", ">=65"), levels = c("<65", ">=65")),
  n = c(107L, 168L, 200L, 150L),
  hr = c(0.61, 0.49, 0.69, 1.19),
  ci_lower = c(0.32, 0.29, 0.45, 0.70),
  ci_upper = c(1.16, 0.84, 1.07, 2.03),
  ci_level = c(0.975, 0.975, 0.951, 0.951),
  stringsAsFactors = FALSE
)
age_tr <- respire_derive(
  respire_age_src$hr,
  respire_age_src$ci_lower,
  respire_age_src$ci_upper,
  respire_age_src$ci_level
)
respire_age <- cbind(respire_age_src, yi = age_tr$yi, sei = age_tr$sei, vi = age_tr$vi)

respire_sex_src <- data.frame(
  trial = c("RESPIRE 1", "RESPIRE 1", "RESPIRE 2", "RESPIRE 2"),
  subgroup = factor(c("Female", "Male", "Female", "Male"), levels = c("Male", "Female")),
  n = c(184L, 91L, 210L, 140L),
  hr = c(0.49, 0.63, 0.64, 1.49),
  ci_lower = c(0.30, 0.31, 0.40, 0.84),
  ci_upper = c(0.81, 1.29, 1.00, 2.63),
  ci_level = c(0.975, 0.975, 0.951, 0.951),
  stringsAsFactors = FALSE
)
sex_tr <- respire_derive(
  respire_sex_src$hr,
  respire_sex_src$ci_lower,
  respire_sex_src$ci_upper,
  respire_sex_src$ci_level
)
respire_sex <- cbind(respire_sex_src, yi = sex_tr$yi, sei = sex_tr$sei, vi = sex_tr$vi)

dir.create("data", showWarnings = FALSE)
save(respire_age, file = "data/respire_age.rda", compress = "xz")
save(respire_sex, file = "data/respire_sex.rda", compress = "xz")
tools::resaveRdaFiles("data", compress = "xz")
