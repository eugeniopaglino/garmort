library(testthat)

test_that("solve_ax_exponential and check_monotonicity behave as expected", {
  # Constant hazard & uniform exposure -> n / 2
  expect_equal(solve_ax_exponential(n = 5, beta = 0, r = 0), 2.5, tolerance = 1e-6)
  # Positive mortality slope -> ax > n / 2
  expect_gt(solve_ax_exponential(n = 5, beta = 0.08, r = 0), 2.5)
  
  # Monotonicity check
  ages <- 0:90
  mx_ok <- 0.0002 * exp(0.08 * ages)
  expect_true(check_monotonicity(ages, mx_ok, min_age = 50)$ok)
  
  mx_bad <- mx_ok
  mx_bad[ages >= 80] <- mx_bad[ages >= 80] * 0.5
  chk <- check_monotonicity(ages, mx_bad, min_age = 50)
  expect_false(chk$ok)
  expect_gt(chk$n_violations, 0L)
})

test_that("abridged-to-single-year reconstruction closely matches true single-year life table", {
  # Construct a realistic Gompertz-Makeham-infant true single-year schedule
  ages_1y <- 0:110
  x_mid_1y <- c(0.2, seq(1.5, 110.5, by = 1))
  mx_true <- 0.00015 + 0.00004 * exp(0.09 * x_mid_1y)
  mx_true[1] <- 0.0045  # infant mortality
  mx_true <- pmin(mx_true, 0.75)
  
  # Build a reference standard (slightly shifted from truth, like HMD vs target country)
  std_raw <- data.table(x = x_mid_1y, mx = mx_true * 1.15)
  std_obj <- build_standard(std_raw, pred_ages = x_mid_1y, kannisto = TRUE)
  expect_true(std_obj$diagnostics$monotone_ok)
  
  # Aggregate true single-year rates into standard abridged counts (0, 1-4, 5-9, ..., 85+)
  abridged_starts <- c(0, 1, seq(5, 85, by = 5))
  grp_1y <- cut(ages_1y, breaks = c(abridged_starts, Inf), right = FALSE, labels = abridged_starts)
  sim_pop <- data.table(age = ages_1y, age_grp = as.numeric(as.character(grp_1y)),
                        Ex = 100000 * exp(-0.015 * ages_1y), mx = mx_true)
  sim_pop[, Dx := round(Ex * mx)]
  
  abridged_counts <- sim_pop[, .(Dx = sum(Dx), Ex = sum(Ex)), by = .(age = age_grp)]
  
  # Fit relational GAM from abridged counts
  fits <- fit_mortality_models(
    data = abridged_counts,
    standard = std_obj$standard,
    age_col = "age", dx_col = "Dx", ex_col = "Ex",
    compute_midpoints = TRUE,
    exclude_open = TRUE,
    pred_ages = x_mid_1y
  )
  
  expect_true(fits$diagnostics$monotone_ok)
  
  preds <- predict_mortality_rates(fits, standard = std_obj$standard, pred_ages = x_mid_1y)
  lt_mod  <- lt_single(data.table(x = ages_1y, mx_s = preds$estimates$mx_s), radix = 1)
  lt_true <- lt_single(data.table(x = ages_1y, mx_s = mx_true), radix = 1)
  
  # "Close match" assertions
  expect_lt(abs(lt_mod[x == 0]$ex  - lt_true[x == 0]$ex),  0.25)
  expect_lt(abs(lt_mod[x == 20]$ex - lt_true[x == 20]$ex), 0.25)
  expect_lt(abs(lt_mod[x == 60]$ex - lt_true[x == 60]$ex), 0.25)
  expect_lt(max(abs(lt_mod[x <= 85]$lx - lt_true[x <= 85]$lx)), 0.01)
  
  log_err <- abs(log(preds$estimates$mx_s[2:86]) - log(mx_true[2:86]))
  expect_lt(median(log_err), 0.05)
})

test_that("sparse/problematic data with 50+ deceleration falls back cleanly", {
  x_mid_1y <- c(0.2, seq(1.5, 110.5, by = 1))
  std_raw <- data.table(x = x_mid_1y, mx = pmin(0.00015 + 0.00004 * exp(0.09 * x_mid_1y), 0.75))
  std_obj <- build_standard(std_raw, pred_ages = x_mid_1y, kannisto = TRUE)
  
  # Sparse data with implausible old-age drop
  abridged_starts <- c(0, 1, seq(5, 85, by = 5))
  prob_dt <- data.table(
    age = abridged_starts,
    Ex  = rep(200, length(abridged_starts))
  )
  prob_dt[, mx := 0.0002 * exp(0.08 * age)]
  prob_dt[age >= 75, mx := mx * 0.3]  # artificial decline
  prob_dt[, Dx := pmin(floor(Ex * mx), Ex - 1)]
  
  fits <- fit_mortality_models(
    data = prob_dt,
    standard = std_obj$standard,
    age_col = "age", dx_col = "Dx", ex_col = "Ex",
    compute_midpoints = TRUE,
    exclude_open = TRUE,
    pred_ages = x_mid_1y
  )
  
  expect_true(fits$diagnostics$monotone_ok)
})