# garmort

<!-- badges: start -->
[![R-CMD-check](https://github.com/eugeniopaglino/garmort/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eugeniopaglino/garmort/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**garmort** provides modular tools for constructing smoothed reference mortality standards, fitting count-based relational generalized additive mortality models (GARM) across demographic strata, propagating parameter uncertainty via empirical Bayes posterior simulation, and computing single-year and abridged life tables.

## Key Features

* **Reference Standard Construction (`build_standard()`):** Smooths reference mortality schedules (e.g., HMD) with penalized thin-plate splines, extrapolates/blends old ages with a Kannisto logistic tail, and enforces non-decreasing mortality past age 50.
* **Count-Based Relational GAM Modeling (`fit_mortality_models()`):** Models death counts (`Dx`) with negative binomial likelihoods and log-exposure offsets (`offset(log(Ex))`) directly across `data.table` grouping columns (`by`).
* **Automated Plausibility Checks & Fallback Ladder:** Evaluates candidate models and verifies convergence plus non-decreasing mortality past age 50 (`monotone_age = 50`). Automatically selects among:
  1. `garm_flexible` — Relational GAM with penalized age spline
  2. `garm_monotonic` — Shape-restricted monotonic spline (`scam`, if installed) or low-`k` penalized GAM
  3. `gam_kannisto` — Relational GAM blended into a Kannisto old-age tail
  4. `garm_joint` — Shared-shape relational GAM pooling sibling strata (e.g., across `nativity`)
  5. `brass` — Two-parameter relational model (`~ 1 + log_mx_std`) for sparse count series
* **Manual Overrides & Per-Group Diagnostics:** Supports explicit per-stratum model overrides (`force_models`) and returns a structured `diagnostics` `data.table` recording selected models, BIC, EDF, total deaths/exposures, and monotonicity check results.
* **Posterior Simulation & Life Tables (`predict_mortality_rates()`, `lt_single()`, `lt_abridged()`, `sim_life_table()`):** Draws from the posterior covariance matrix (`Vp`) to produce single-year and abridged life tables with 95% credible intervals, plus analytical exponential-hazard $a_x$ solvers (`solve_ax_exponential()`).
* **Cross-Platform Parallelization:** Supports multi-worker execution on macOS, Windows, and Linux via `future.apply`.

## Installation

Install the development version directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("<your-username>/garmort")
```

Optional packages enabling monotonic spline fallbacks, cross-platform parallel workers, and diagnostic plots:

```r
install.packages(c("scam", "future", "future.apply", "ggplot2"))
```

## Quick Start

```r
library(data.table)
library(garmort)

# 1. Build a smoothed reference standard (e.g., from HMD rates)
std_res <- build_standard(
  data = ref_rates,       # columns: sex, x, mx
  age_col = "x",
  mx_col = "mx",
  by = "sex",
  kannisto = TRUE
)

# 2. Fit relational mortality models across groups using abridged counts
fits <- fit_mortality_models(
  data = observed_counts, # columns: state, year, sex, nativity, age, Dx, Ex
  standard = std_res$standard,
  by = c("state", "year", "sex", "nativity"),
  std_by = "sex",
  age_col = "age",
  dx_col = "Dx",
  ex_col = "Ex",
  compute_midpoints = TRUE,
  exclude_open = TRUE,
  joint_col = "nativity",
  workers = 4L
)

# Inspect model selection and 50+ monotonicity diagnostics per series
print(fits$diagnostics)

# 3. Predict single-year rates and draw posterior simulations
preds <- predict_mortality_rates(
  object = fits,
  standard = std_res$standard,
  pred_ages = c(0.2, seq(1.5, 110.5, by = 1)),
  n_sim = 500L,
  seed = 123L
)

# 4. Compute point-estimate or posterior-summarized life tables
est <- copy(preds$estimates)[, age := fifelse(x == 0.2, 0L, as.integer(floor(x)))]
lt_1y <- lt_single(est, age_col = "age", mx_col = "mx_s",
                   by = c("state", "year", "sex", "nativity"))

sim_lts <- sim_life_table(
  simulations = preds$simulations,
  by = c("state", "year", "sex", "nativity"),
  max_abridged_age = 85,
  radix = 1
)
```

## Core Functions

| Function | Purpose |
| :--- | :--- |
| `build_standard()` | Smooth reference mortality rates and blend Kannisto old-age tails with 50+ monotonicity checks |
| `fit_mortality_models()` | Fit grouped negative binomial relational GAMs with automated fallback ladder and diagnostics |
| `predict_mortality_rates()` | Predict point-estimate rates (`mx_s`) and posterior rate draws (`sim_1..sim_n`) |
| `lt_single()` / `lt_abridged()` | Compute single-year and abridged life tables across grouped `data.table`s |
| `sim_life_table()` | Aggregate posterior single-year rate simulations into single-year and abridged life table summaries with credible intervals |
| `compute_age_intervals()` | Map variable starting/ending age schedules to interval widths and midpoints |
| `solve_ax_exponential()` | Analytical $a_x$ solver under exponential mortality and exposure slopes |
| `plot_mortality_fit()` | Visual diagnostic comparing observed, modeled, and standard log-mortality schedules |
