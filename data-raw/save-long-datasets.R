# Write long-format package datasets (from current wide .rda on disk).
# Run from repo root: Rscript data-raw/save-long-datasets.R

root <- normalizePath(".", wins = FALSE)
source(file.path(root, "R", "subgroup_long.R"), local = TRUE)

for (nm in c("age", "edss", "il6", "steroids", "respire_age", "respire_sex")) {
  load(file.path(root, "data", paste0(nm, ".rda")), envir = .GlobalEnv)
  d <- get(nm)
  long <- if ("subgroup" %in% names(d)) d else subgroup_long(d)
  long$subgroup <- as.character(long$subgroup)
  long <- switch(
    nm,
    age = long[, c("study", "subgroup", "n", "treatment", "yi", "vi")],
    edss = long[, c("study", "subgroup", "n", "treatment", "yi", "vi")],
    il6 = long[, c("study", "subgroup", "n", "treat.events", "treat.total", "cont.events", "cont.total")],
    steroids = long[, c("study", "subgroup", "treatment", "n", "treat.events", "treat.total", "cont.events", "cont.total")],
    respire_age = long,
    respire_sex = long,
    long
  )
  if (nm == "edss") {
    long$subgroup <- factor(long$subgroup, levels = c("Lower", "Higher"))
  }
  if (nm == "respire_age") {
    long$subgroup <- factor(long$subgroup, levels = c("<65", ">=65"))
  }
  if (nm == "respire_sex") {
    long$subgroup <- factor(long$subgroup, levels = c("Female", "Male"))
  }
  assign(nm, long, envir = .GlobalEnv)
  save(list = nm, file = file.path(root, "data", paste0(nm, ".rda")), compress = "xz")
  message("Wrote long ", nm, " (", nrow(long), " rows)")
}
