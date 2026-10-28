# interact test tiers

## Fast (default)

Used by `R CMD check`, PR CI, and local dev:

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Slow plot and full forestplot integration tests are **skipped** unless you opt in.

## Slow tier

Set `INTERACT_SLOW_TESTS=true` for end-to-end forestplot PDF/render coverage:

```bash
INTERACT_SLOW_TESTS=true Rscript -e 'testthat::test_dir("tests/testthat")'
```

Files that call `slow()` at the top of a test include:

- [`test-forestplot-integration.R`](test-forestplot-integration.R)
- [`test-forestplot_comparison.R`](test-forestplot_comparison.R) (render paths)
- [`test-shrinkage_forest.R`](test-shrinkage_forest.R)

## Fixtures

Shared data and cached fits live in [`helper-fixtures.R`](helper-fixtures.R).

## Coverage and timing

```bash
Rscript scripts/test-timing.R
INTERACT_SLOW_TESTS=true Rscript scripts/coverage.R
```

HTML coverage is written to `coverage/` (gitignored).
