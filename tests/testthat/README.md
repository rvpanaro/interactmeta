# interactmeta tests

## Fast (default)

Used by `R CMD check`, PR CI, and local dev.

**RStudio / Console:**

```r
devtools::test()
```

Or:

```r
pkgload::load_all()
testthat::test_dir("tests/testthat")
```

Helpers in `helper-fixtures.R` call `pkgload::load_all()` when the package is not attached.

```bash
Rscript -e 'pkgload::load_all(); testthat::test_dir("tests/testthat")'
```

## Suite layout

| File | Coverage |
|------|----------|
| `test-battery-cams.R` | `cams()` backends, removed args, lazy precompute |
| `test-battery-ifrac.R` | Design tables and `ifrac()` on data / fits |
| `test-battery-effects.R` | `modifiers`, `interactions`, `combined`, `heterogeneity`, `shrinkage` |
| `test-battery-constructions.R` | `da`, `ad`, `swada`, `conventional` (REML) |
| `test-battery-reporting.R` | `print`, `summary`, `tidy` |
| `test-battery-compare.R` | `compare`, `compare_weights`, `sensitivity()`, capability matrix |
| `test-battery-data.R` | Datasets and `log_hr_from_ci` / `hr_ci_from_log` |
| `test-battery-graphics.R` | `forestplot`, `forestplot_draw`, `traceplot`, `compare_forestplots`, plots (`slow()` tier) |

## Slow tier

Several graphics / integration tests are skipped unless you opt in (`slow()` helper):

```bash
INTERACT_SLOW_TESTS=true Rscript -e 'pkgload::load_all(); testthat::test_dir("tests/testthat")'
```

## Fixtures

Shared data and cached fits: `helper-fixtures.R`.

## Coverage

```bash
INTERACT_SLOW_TESTS=true Rscript scripts/coverage.R
```

HTML coverage is written to `coverage/` (gitignored).