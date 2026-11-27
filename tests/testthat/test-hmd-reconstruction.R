library(testthat)
library(data.table)

test_that("HMD abridged counts reconstruct single-year rates and life tables accurately", {
  rds_path <- system.file("extdata", "hmd_test_data.rds", package = "garmort")
  skip_if(rds_path == "", message = "HMD test fixture not found in inst/extdata.")
  
  hmd <- readRDS(rds_path)
  
  # Include x = 3 so the [1, 5) abridged midpoint matches during model fitting
  std_grid <- sort(unique(c(0.2, 3, seq(1.5, 110.5, by = 1))))
  pred_1y  <- c(0.2, seq(1.5, 110.5, by = 1))
  
  std_res <- build_standard(
    data = hmd$std_input,
    age_col = "x",
    mx_col = "mx",
    by = "sex",
    pred_ages = std_grid,
    kannisto = TRUE
  )
  expect_true(all(std_res$diagnostics$monotone_ok))
  
  # Fit across all 4 country-years x 3 sexes (12 series) simultaneously
  fits <- fit_mortality_models(
    data = hmd$abridged,
    standard = std_res$standard,
    by = c("country", "year", "sex"),
    std_by = "sex",
    age_col = "age",
    dx_col = "Dx",
    ex_col = "Ex",
    compute_midpoints = TRUE,
    exclude_open = TRUE,
    pred_ages = pred_1y
  )
  
  # Every reconstructed curve must obey 50+ monotonicity
  expect_true(all(fits$diagnostics$monotone_ok))
  
  # Point predictions + posterior simulations
  preds <- predict_mortality_rates(
    fits,
    standard = std_res$standard,
    pred_ages = pred_1y,
    n_sim = 300L,
    seed = 42L
  )
  
  est <- copy(preds$estimates)
  est[, age := fifelse(x == 0.2, 0L, as.integer(floor(x)))]
  
  # Point-estimate single-year life tables
  lt_est <- lt_single(est, age_col = "age", mx_col = "mx_s",
                      by = c("country", "year", "sex"), radix = 1)
  
  # Posterior simulated single-year & abridged life table summaries
  sim_lt <- sim_life_table(
    preds$simulations,
    by = c("country", "year", "sex"),
    max_abridged_age = 85,
    radix = 1
  )
  s_sum <- sim_lt$single_year_summary
  
  # Verify credible interval ordering
  expect_true(all(s_sum$mx_l95 <= s_sum$mx_mean & s_sum$mx_mean <= s_sum$mx_u95))
  expect_true(all(s_sum$ex_l95 <= s_sum$ex_mean & s_sum$ex_mean <= s_sum$ex_u95))
  
  comp <- merge(
    lt_est[, .(country, year, sex, age, mx_s = mx, lx_s = lx, ex_s = ex)],
    hmd$truth_1y,
    by = c("country", "year", "sex", "age")
  )
  
  # Accuracy checks by country-year (Total sex)
  tot <- comp[sex == "Total"]
  
  # 1. Life expectancy at birth (e0) and age 60 (e60)
  e0_err <- tot[age == 0, .(err = abs(ex_s - ex_true)), by = country]
  expect_lt(e0_err[country == "USA"]$err, 0.35)
  expect_lt(e0_err[country == "RUS"]$err, 0.35)
  expect_lt(e0_err[country == "JPN"]$err, 0.35)
  expect_lt(e0_err[country == "ISL"]$err, 0.60)
  
  e60_err <- tot[age == 60, .(err = abs(ex_s - ex_true)), by = country]
  expect_lt(max(e60_err$err), 0.35)
  
  # 2. Single-year log-rate median absolute error across ages 1-84
  rate_err <- tot[age %between% c(1L, 84L) & mx_true > 0,
                  .(med_log_err = median(abs(log(mx_s) - log(mx_true)))),
                  by = country]
  expect_lt(rate_err[country == "USA"]$med_log_err, 0.05)
  expect_lt(rate_err[country == "RUS"]$med_log_err, 0.06)
  expect_lt(rate_err[country == "JPN"]$med_log_err, 0.06)
  expect_lt(rate_err[country == "ISL"]$med_log_err, 0.25)
  
  # 3. Save diagnostic PDFs with 95% credible interval ribbons
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_dir <- get_test_plot_dir()
    
    obs_plot <- copy(hmd$abridged)
    obs_plot[, c("x", "is_open") := {
      iv <- compute_age_intervals(age)
      idx <- match(age, iv$age)
      list(iv$x_mid[idx], iv$is_open[idx])
    }, by = .(country, year, sex)]
    obs_plot[, mx_obs := Dx / Ex]
    
    diag_labels <- fits$diagnostics[, .(
      country, year, sex,
      label = paste0("Model: ", selected_model, " (deaths=", round(total_deaths), ")")
    )]
    
    comp_plot <- copy(comp)
    comp_plot[, x := fifelse(age == 0L, 0.2, age + 0.5)]
    
    s_sum_plot <- copy(s_sum)
    s_sum_plot[, x_mid := fifelse(x == 0L, 0.2, x + 0.5)]
    
    # Plot A: Log-mortality rates + 95% CI ribbon
    p_mx <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = s_sum_plot,
        ggplot2::aes(x = x_mid, ymin = log(mx_l95), ymax = log(mx_u95), fill = "95% Credible Interval"),
        alpha = 0.35
      ) +
      ggplot2::geom_line(
        data = preds$estimates,
        ggplot2::aes(x = x, y = log(mx_std), color = "Standard"),
        linetype = "dotted", linewidth = 0.8
      ) +
      ggplot2::geom_line(
        data = comp_plot[mx_true > 0],
        ggplot2::aes(x = x, y = log(mx_true), color = "True 1x1 HMD"),
        linetype = "dashed", linewidth = 0.8
      ) +
      ggplot2::geom_line(
        data = preds$estimates,
        ggplot2::aes(x = x, y = log(mx_s), color = "Modeled 1x1"),
        linewidth = 0.9
      ) +
      ggplot2::geom_point(
        data = obs_plot[mx_obs > 0],
        ggplot2::aes(x = x, y = log(mx_obs), shape = is_open),
        color = "black", size = 1.6, alpha = 0.8
      ) +
      ggplot2::geom_text(
        data = diag_labels,
        ggplot2::aes(x = 5, y = -0.5, label = label),
        hjust = 0, vjust = 1, size = 3, color = "grey20"
      ) +
      ggplot2::facet_grid(country ~ sex) +
      ggplot2::scale_fill_manual(values = c("95% Credible Interval" = "#E29595")) +
      ggplot2::scale_shape_manual(
        values = c("FALSE" = 16, "TRUE" = 1),
        labels = c("FALSE" = "Fitted Abridged", "TRUE" = "Excluded Open (85+)")
      ) +
      ggplot2::scale_color_manual(
        values = c("Modeled 1x1" = "#B22222", "True 1x1 HMD" = "black", "Standard" = "#2B5C8F")
      ) +
      ggplot2::labs(
        title = "HMD Abridged-to-Single-Year Reconstruction: Log Mortality Rates",
        x = "Age", y = "log(mx)", color = NULL, fill = NULL, shape = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    
    grDevices::pdf(file.path(plot_dir, "hmd_log_mx_fit.pdf"), width = 12, height = 10)
    print(p_mx)
    grDevices::dev.off()
    
    # Plot B: Remaining life expectancy (ex) + 95% CI ribbon
    e0_labels <- comp[age == 0, .(
      label = sprintf("e0 diff: %+.2f yrs", ex_s - ex_true)
    ), by = .(country, year, sex)]
    
    p_ex <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = s_sum[x <= 100],
        ggplot2::aes(x = x, ymin = ex_l95, ymax = ex_u95, fill = "95% Credible Interval"),
        alpha = 0.35
      ) +
      ggplot2::geom_line(
        data = comp[age <= 100],
        ggplot2::aes(x = age, y = ex_true, color = "True 1x1 HMD"),
        linetype = "dashed", linewidth = 0.9
      ) +
      ggplot2::geom_line(
        data = comp[age <= 100],
        ggplot2::aes(x = age, y = ex_s, color = "Modeled 1x1"),
        linewidth = 0.9
      ) +
      ggplot2::geom_text(
        data = e0_labels,
        ggplot2::aes(x = 58, y = 78, label = label),
        hjust = 0, size = 3.2, color = "grey20"
      ) +
      ggplot2::facet_grid(country ~ sex) +
      ggplot2::scale_fill_manual(values = c("95% Credible Interval" = "#E29595")) +
      ggplot2::scale_color_manual(
        values = c("Modeled 1x1" = "#B22222", "True 1x1 HMD" = "black")
      ) +
      ggplot2::labs(
        title = "HMD Abridged-to-Single-Year Reconstruction: Remaining Life Expectancy (ex)",
        x = "Exact Age", y = "Life Expectancy (Years)", color = NULL, fill = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    
    grDevices::pdf(file.path(plot_dir, "hmd_ex_comparison.pdf"), width = 12, height = 10)
    print(p_ex)
    grDevices::dev.off()
  }
})