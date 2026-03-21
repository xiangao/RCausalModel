# RCausalModel

## Overview
R package port of Python causal inference library implementing estimators for observational studies, randomized experiments, and network interference settings. Based on Qu et al. (2021). Python is the reference implementation.

## Structure
- `R/` — package source (interference.R, random_data.R, observational.R, etc.)
- `tests/testthat/` — 113 tests including slow simulation tests
- `vignettes/` — `quickstart.Rmd` (usage walkthrough), `simulations.Rmd` (Monte Carlo validation), `observational.Rmd` (mirrors Python observational.ipynb), `experimental.Rmd` (mirrors Python experimental.ipynb), `interference.Rmd` (mirrors Python Interference.ipynb)
- Parent project with Python source: `~/projects/claude/CausalModel`

## Key Commands
```bash
Rscript -e 'devtools::test()'
Rscript -e 'devtools::check(cran=FALSE)'
Rscript -e 'devtools::document()'
Rscript -e 'devtools::build_vignettes()'
```

## Testing
- Some simulation tests (normality, coverage) are slow (10+ minutes) — they are NOT skipped
- `skip_on_cran()` is used for slow tests but they run locally
- All estimators: OLS, IPW, AIPW, matching, DML, difference-in-means, strata, ANCOVA, Fisher test, clustered IPW/AIPW

## Design Decisions
- S3 classes (idiomatic R) rather than R6/R5
- Model wrappers use closure-based protocol: list with `fit(X, y)` and `predict(model, X)` functions
- HC0 robust SEs computed manually to avoid `sandwich` dependency
- RANN::nn2() replaces scipy.spatial.KDTree for matching
- Dependencies: stats (base), RANN (matching), nnet (multinomial logistic)

## Gotchas
- R is column-major: when flattening 3D arrays (clusters, nunit, ngroup) to 2D matrices, use `aperm(arr, c(2, 1, 3))` before `matrix()` to match Python's row-major reshape behavior
- Python is the reference implementation — when R results diverge, the R code is likely wrong
- Watch for R/Python differences: column-major vs row-major, 1-indexing vs 0-indexing

## Bug Fixes
- **interference.R `get_final_tuple` aperm fix** (2026-03-17): The augmented covariate matrix X_aug (3D: clusters × nunit × k_aug) was flattened to 2D via `matrix(X_aug, ...)` without `aperm()`, scrambling covariate rows relative to Y, Z, G. Fix: `matrix(aperm(X_aug, c(2,1,3)), ...)`. This was the same column-major reshape issue documented in Gotchas — the G computation already had the fix but X did not. Impact: AIPW coverage at rare g values went from ~60% to ~94%.
