#' Build a Smoothed Reference Mortality Standard
#'
#' Smooths reference mortality rates (e.g., from HMD) over age by group,
#' optionally applies Kannisto old-age extrapolation and blending, and verifies
#' that mortality is non-decreasing past `monotone_age`.
#'
#' @param data A `data.table` (or data.frame) containing reference mortality rates.
#' @param age_col Name of the numeric age column (default `"x"`).
#' @param mx_col Name of the mortality rate column (default `"mx"`).
#' @param by Optional character vector of grouping columns (e.g., `"sex"` or `c("year", "sex")`).
#' @param pred_ages Numeric vector of ages to predict on (default `sort(unique(c(0.2, 3, seq(1.5, 111.5, by = 1))))`).
#' @param k Integer spline basis dimension for the GAM (default 30).
#' @param kannisto Logical. Apply Kannisto extrapolation and blending at high ages? (default TRUE).
#' @param kannisto_fit_ages Numeric length-2 vector for Kannisto fit window (default `c(80, 105)`).
#' @param kannisto_blend_ages Numeric length-2 vector for linear blending window (default `c(100, 105)`).
#' @param monotone_age Numeric age above which predicted rates must be non-decreasing (default 50).
#' @return A list with:
#'   - `standard`: `data.table` with grouping columns, `x`, `mx_std`, `log_mx_std`.
#'   - `diagnostics`: `data.table` with per-group method used and monotonicity check results.
#' @importFrom data.table as.data.table copy rbindlist setnames fcase
#' @importFrom mgcv gam s
#' @importFrom stats lm predict qlogis plogis
#' @export
build_standard <- function(data,
                           age_col = "x",
                           mx_col = "mx",
                           by = NULL,
                           pred_ages = sort(unique(c(0.2, 3, seq(1.5, 111.5, by = 1)))),
                           k = 30L,
                           kannisto = TRUE,
                           kannisto_fit_ages = c(80, 105),
                           kannisto_blend_ages = c(100, 105),
                           monotone_age = 50) {
  dt <- data.table::copy(data.table::as.data.table(data))
  if (!age_col %in% names(dt)) stop("Column `", age_col, "` not found in `data`.")
  if (!mx_col %in% names(dt)) stop("Column `", mx_col, "` not found in `data`.")
  
  dt[, .x := as.numeric(get(age_col))]
  dt[, .log_mx := log(as.numeric(get(mx_col)))]
  dt <- dt[is.finite(.x) & is.finite(.log_mx)]
  
  if (is.null(by) || length(by) == 0L) {
    dt[, .grp_id := 1L]
    grp_keys <- data.table::data.table(.grp_id = 1L)
  } else {
    dt[, .grp_id := .GRP, by = by]
    grp_keys <- unique(dt[, c(".grp_id", by), with = FALSE])
  }
  
  fit_one_std <- function(sub_dt) {
    k_use <- min(as.integer(k), length(unique(sub_dt$.x)) - 2L)
    method_used <- "gam_tp"
    
    fit <- mgcv::gam(.log_mx ~ s(.x, k = k_use, bs = "tp"),
                     data = sub_dt, method = "REML")
    pred_dt <- data.table::data.table(.x = pred_ages)
    pred_dt[, log_mx_smooth := as.numeric(stats::predict(fit, newdata = pred_dt))]
    
    # Check monotonicity past monotone_age on smoothed curve
    chk <- check_monotonicity(pred_dt$.x, pred_dt$log_mx_smooth, min_age = monotone_age)
    
    if (!chk$ok) {
      # Fallback 1: lower-k penalized GAM
      k_low <- max(6L, as.integer(k_use / 3))
      fit_low <- mgcv::gam(.log_mx ~ s(.x, k = k_low, bs = "tp"),
                           data = sub_dt, method = "REML")
      pred_low <- as.numeric(stats::predict(fit_low, newdata = pred_dt))
      chk_low <- check_monotonicity(pred_dt$.x, pred_low, min_age = monotone_age)
      
      if (chk_low$ok) {
        pred_dt[, log_mx_smooth := pred_low]
        method_used <- "gam_tp_low_k"
        chk <- chk_low
      } else if (requireNamespace("scam", quietly = TRUE)) {
        # Fallback 2: monotonic increasing tail via scam on >= monotone_age
        fit_scam <- tryCatch(
          scam::scam(.log_mx ~ s(.x, k = k_low, bs = "mpi"), data = sub_dt[.x >= monotone_age]),
          error = function(e) NULL
        )
        if (!is.null(fit_scam)) {
          idx_hi <- pred_dt$.x >= monotone_age
          pred_dt[idx_hi, log_mx_smooth := as.numeric(stats::predict(fit_scam, newdata = pred_dt[idx_hi]))]
          method_used <- "gam_plus_scam_mpi"
          chk <- check_monotonicity(pred_dt$.x, pred_dt$log_mx_smooth, min_age = monotone_age)
        }
      }
    }
    
    pred_dt[, mx_smooth := exp(log_mx_smooth)]
    
    # Optional Kannisto extrapolation & blending
    if (isTRUE(kannisto)) {
      fit_k <- pred_dt[.x >= kannisto_fit_ages[1] & .x <= kannisto_fit_ages[2]]
      if (nrow(fit_k) >= 5L) {
        fit_k[, mx_cap := pmin(pmax(mx_smooth, 1e-6), 0.999)]
        mod_k <- stats::lm(stats::qlogis(mx_cap) ~ I(.x - kannisto_fit_ages[1]), data = fit_k)
        pred_dt[.x >= kannisto_fit_ages[1],
                mx_kannisto := stats::plogis(stats::predict(mod_k, newdata = .SD))]
        
        b1 <- kannisto_blend_ages[1]
        b2 <- kannisto_blend_ages[2]
        pred_dt[, w_k := data.table::fcase(
          .x < b1, 0,
          .x >= b1 & .x < b2, (.x - b1) / (b2 - b1),
          default = 1
        )]
        pred_dt[, mx_std := data.table::fcase(
          w_k == 0, mx_smooth,
          w_k == 1, mx_kannisto,
          default = (1 - w_k) * mx_smooth + w_k * mx_kannisto
        )]
      } else {
        pred_dt[, mx_std := mx_smooth]
      }
    } else {
      pred_dt[, mx_std := mx_smooth]
    }
    
    pred_dt[, log_mx_std := log(mx_std)]
    final_chk <- check_monotonicity(pred_dt$.x, pred_dt$log_mx_std, min_age = monotone_age)
    
    list(
      std = pred_dt[, .(x = .x, mx_std, log_mx_std)],
      diag = data.table::data.table(
        method = method_used,
        kannisto = isTRUE(kannisto),
        monotone_ok = final_chk$ok,
        n_violations = final_chk$n_violations,
        max_drop = final_chk$max_drop
      )
    )
  }
  
  std_list <- vector("list", nrow(grp_keys))
  diag_list <- vector("list", nrow(grp_keys))
  
  for (i in seq_len(nrow(grp_keys))) {
    gid <- grp_keys$.grp_id[i]
    res <- fit_one_std(dt[.grp_id == gid])
    if (!is.null(by) && length(by) > 0L) {
      for (col in by) {
        res$std[, (col) := grp_keys[[col]][i]]
        res$diag[, (col) := grp_keys[[col]][i]]
      }
    }
    std_list[[i]] <- res$std
    diag_list[[i]] <- res$diag
  }
  
  std_out <- data.table::rbindlist(std_list, use.names = TRUE)
  diag_out <- data.table::rbindlist(diag_list, use.names = TRUE)
  
  if (!is.null(by) && length(by) > 0L) {
    data.table::setcolorder(std_out, c(by, "x", "mx_std", "log_mx_std"))
    data.table::setcolorder(diag_out, c(by, "method", "kannisto", "monotone_ok", "n_violations", "max_drop"))
  }
  
  list(standard = std_out, diagnostics = diag_out)
}