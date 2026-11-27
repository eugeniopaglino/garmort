#' Fit Count-Based Relational GAM Mortality Models Across Groups
#'
#' Fits negative binomial relational mortality models to death counts (`Dx`) and
#' exposures (`Ex`) against a reference standard, with automated BIC selection,
#' post-fit demographic plausibility checks (non-decreasing mortality past
#' `monotone_age`), optional joint shared-shape pooling across sibling strata,
#' manual model overrides, and cross-platform parallel execution.
#'
#' @param data A `data.table` (or data.frame) of observed mortality counts.
#' @param standard A `data.table` (or data.frame) of reference standard rates
#'   (e.g. from `build_standard()$standard`), containing `x` and `log_mx_std`
#'   (or `mx_std`).
#' @param by Character vector of grouping columns defining each series
#'   (e.g., `c("state_abb", "year", "sex", "nativity")`).
#' @param std_by Character vector of columns used to join `standard` to `data`
#'   (default `intersect(by, names(standard))`).
#' @param age_col Name of the age column in `data` (default `"x"`). Can be
#'   interval start ages or midpoints; if `compute_midpoints = TRUE`, midpoints
#'   are computed automatically via `compute_age_intervals()`.
#' @param dx_col Name of the death count column (default `"Dx"`).
#' @param ex_col Name of the exposure/population column (default `"Ex"`).
#' @param compute_midpoints Logical. If `TRUE`, treats `age_col` as interval
#'   starting ages and maps them to midpoints before joining `standard` (default `FALSE`).
#' @param exclude_open Logical. Exclude the highest (open-ended) age group from
#'   model fitting? (default `TRUE`).
#' @param skip_ages Optional list of length-2 numeric vectors `c(low, high)`
#'   defining additional age ranges to exclude from fitting.
#' @param pred_ages Numeric vector of ages used for post-fit plausibility checks
#'   and default prediction (default `c(0.2, seq(1.5, 111.5, by = 1))`).
#' @param dof Integer spline basis dimension for flexible GARM (default `10L`).
#' @param basis Spline basis type for `mgcv::s()` (default `"tp"`).
#' @param min_deaths_for_smooth Minimum total deaths required to attempt a smooth
#'   age term before falling back to Brass (default `250`).
#' @param monotone_age Numeric age above which predicted mortality must be
#'   non-decreasing (default `50`).
#' @param joint_col Optional single column name in `by` across which sibling
#'   strata can share a common smooth age shape (e.g., `"nativity"`).
#' @param joint_exclude Optional character vector of values in `joint_col` to
#'   exclude from joint pooling (e.g., `"total"` or `"Total"`).
#' @param force_models Optional `data.table` specifying manual overrides. Must
#'   contain a `forced_model` column (one of `"garm_flexible"`, `"garm_monotonic"`,
#'   `"gam_kannisto"`, `"garm_joint"`, `"brass"`) and one or more columns from `by`.
#' @param workers Integer number of parallel workers (default `1L`). Values > 1
#'   use `future.apply` with a cross-platform `multisession` plan.
#' @return An S3 object of class `"garmort_fits"` containing:
#'   - `fits`: Named list of fitted model objects (one per group).
#'   - `diagnostics`: `data.table` with per-group model selection and check results.
#'   - `group_dict`: `data.table` mapping `.grp_id` to `by` columns.
#'   - `meta`: List of metadata (`by`, `std_by`, `pred_ages`, `monotone_age`, etc.).
#' @importFrom data.table as.data.table copy rbindlist setnames setorder fcase
#' @importFrom mgcv gam s nb
#' @importFrom stats BIC coef predict lm qlogis plogis quasipoisson
#' @export
fit_mortality_models <- function(data,
                                 standard,
                                 by = NULL,
                                 std_by = NULL,
                                 age_col = "x",
                                 dx_col = "Dx",
                                 ex_col = "Ex",
                                 compute_midpoints = FALSE,
                                 exclude_open = TRUE,
                                 skip_ages = NULL,
                                 pred_ages = c(0.2, seq(1.5, 111.5, by = 1)),
                                 dof = 10L,
                                 basis = "tp",
                                 min_deaths_for_smooth = 250,
                                 monotone_age = 50,
                                 joint_col = NULL,
                                 joint_exclude = c("Total", "total"),
                                 force_models = NULL,
                                 workers = 1L) {
  
  dt <- data.table::copy(data.table::as.data.table(data))
  std <- data.table::copy(data.table::as.data.table(standard))
  
  # Validate required columns
  for (col in c(age_col, dx_col, ex_col)) {
    if (!col %in% names(dt)) stop("Column `", col, "` not found in `data`.")
  }
  if (!"x" %in% names(std)) stop("`standard` must contain column `x`.")
  if (!"log_mx_std" %in% names(std)) {
    if ("mx_std" %in% names(std)) {
      std[, log_mx_std := log(mx_std)]
    } else {
      stop("`standard` must contain `log_mx_std` or `mx_std`.")
    }
  }
  
  # Standardize working column names
  dt[, .age_raw := as.numeric(get(age_col))]
  dt[, Dx := as.numeric(get(dx_col))]
  dt[, Ex := as.numeric(get(ex_col))]
  
  # Grouping setup
  if (is.null(by) || length(by) == 0L) {
    dt[, .grp_id := 1L]
    group_dict <- data.table::data.table(.grp_id = 1L)
  } else {
    dt[, .grp_id := .GRP, by = by]
    group_dict <- unique(dt[, c(".grp_id", by), with = FALSE])
  }
  
  # Optional midpoint computation per group
  if (isTRUE(compute_midpoints)) {
    dt[, c(".x", ".is_open") := {
      iv <- compute_age_intervals(.age_raw)
      idx <- match(.age_raw, iv$age)
      list(iv$x_mid[idx], iv$is_open[idx])
    }, by = .grp_id]
  } else {
    dt[, .x := .age_raw]
    dt[, .is_open := (.x == max(.x, na.rm = TRUE)), by = .grp_id]
  }
  
  # Join standard
  if (is.null(std_by)) {
    std_by <- intersect(by, names(std))
  }
  join_cols <- c(std_by, "x")
  dt[, x := .x]
  dt <- merge(dt, std[, c(join_cols, "log_mx_std"), with = FALSE],
              by = join_cols, all.x = TRUE)
  
  # Apply fitting exclusions
  dt[, .use_fit := is.finite(Dx) & is.finite(Ex) & Ex > 0 & Ex > Dx & is.finite(log_mx_std)]
  if (isTRUE(exclude_open)) {
    dt[.is_open == TRUE, .use_fit := FALSE]
  }
  if (!is.null(skip_ages)) {
    if (!is.list(skip_ages)) skip_ages <- list(skip_ages)
    for (rng in skip_ages) {
      if (is.numeric(rng) && length(rng) == 2L) {
        dt[x >= rng[1] & x <= rng[2], .use_fit := FALSE]
      }
    }
  }
  
  # Resolve manual overrides onto group_dict
  group_dict[, forced_model := NA_character_]
  if (!is.null(force_models)) {
    fm <- data.table::as.data.table(force_models)
    if (!"forced_model" %in% names(fm)) {
      stop("`force_models` must contain a `forced_model` column.")
    }
    ov_keys <- intersect(names(fm), by)
    if (length(ov_keys) == 0L) {
      stop("`force_models` must share at least one column with `by`.")
    }
    group_dict[fm, forced_model := i.forced_model, on = ov_keys]
  }
  
  # Pre-build joint cluster IDs if joint_col is requested
  has_joint <- !is.null(joint_col) && length(joint_col) == 1L && (joint_col %in% by)
  if (has_joint) {
    parent_by <- setdiff(by, joint_col)
    if (length(parent_by) > 0L) {
      group_dict[, .joint_id := .GRP, by = parent_by]
    } else {
      group_dict[, .joint_id := 1L]
    }
  } else {
    group_dict[, .joint_id := .grp_id]
  }
  
  # Helper: evaluate candidate model & check 50+ monotonicity on pred_grid
  eval_candidate <- function(fit_obj, pred_grid, model_name, kannisto_tail = FALSE) {
    if (is.null(fit_obj) || inherits(fit_obj, "error")) {
      return(list(ok = FALSE, converged = FALSE, monotone_ok = FALSE,
                  n_violations = NA_integer_, bic = NA_real_, edf = NA_real_,
                  fit = NULL, model = model_name))
    }
    
    conv <- if (!is.null(fit_obj$converged)) isTRUE(fit_obj$converged) else TRUE
    preds <- tryCatch({
      as.numeric(stats::predict(fit_obj, newdata = pred_grid, type = "link"))
    }, error = function(e) NULL)
    
    if (is.null(preds) || any(!is.finite(preds))) {
      return(list(ok = FALSE, converged = FALSE, monotone_ok = FALSE,
                  n_violations = NA_integer_, bic = NA_real_, edf = NA_real_,
                  fit = NULL, model = model_name))
    }
    
    # Optional Kannisto tail blending for gam_kannisto
    kannisto_mod <- NULL
    if (isTRUE(kannisto_tail)) {
      mx_p <- exp(preds)
      k_idx <- which(pred_grid$x >= 80 & pred_grid$x <= 105)
      if (length(k_idx) >= 5L) {
        k_dt <- data.table::data.table(
          x = pred_grid$x[k_idx],
          mx_cap = pmin(pmax(mx_p[k_idx], 1e-6), 0.999)
        )
        kannisto_mod <- tryCatch(
          stats::lm(stats::qlogis(mx_cap) ~ I(x - 80), data = k_dt),
          error = function(e) NULL
        )
        if (!is.null(kannisto_mod)) {
          hi_idx <- which(pred_grid$x >= 80)
          mx_k <- stats::plogis(stats::predict(kannisto_mod, newdata = pred_grid[hi_idx]))
          w_k <- data.table::fcase(
            pred_grid$x[hi_idx] < 100, 0,
            pred_grid$x[hi_idx] >= 100 & pred_grid$x[hi_idx] < 105, (pred_grid$x[hi_idx] - 100) / 5,
            default = 1
          )
          mx_blend <- (1 - w_k) * mx_p[hi_idx] + w_k * mx_k
          preds[hi_idx] <- log(mx_blend)
        }
      }
    }
    
    chk <- check_monotonicity(pred_grid$x, preds, min_age = monotone_age)
    bic_val <- tryCatch(stats::BIC(fit_obj), error = function(e) NA_real_)
    edf_val <- if (!is.null(fit_obj$edf)) sum(fit_obj$edf) else NA_real_
    
    wrapped_fit <- list(
      engine = fit_obj,
      model = model_name,
      kannisto_mod = kannisto_mod
    )
    
    list(
      ok = conv && chk$ok,
      converged = conv,
      monotone_ok = chk$ok,
      n_violations = chk$n_violations,
      bic = bic_val,
      edf = edf_val,
      fit = wrapped_fit,
      model = model_name
    )
  }
  
  # Core per-group fitting function
  fit_one_group <- function(gid) {
    g_row <- group_dict[.grp_id == gid]
    sub_dt <- dt[.grp_id == gid & .use_fit == TRUE]
    total_deaths <- sum(sub_dt$Dx, na.rm = TRUE)
    total_exp <- sum(sub_dt$Ex, na.rm = TRUE)
    
    # Build prediction grid for this group's standard
    if (length(std_by) > 0L) {
      std_sub <- merge(
        data.table::data.table(x = pred_ages),
        std[g_row[, std_by, with = FALSE], on = std_by, nomatch = NULL],
        by = "x"
      )
    } else {
      std_sub <- std[x %in% pred_ages]
    }
    std_sub[, Ex := 1]
    if (has_joint) {
      std_sub[, (joint_col) := g_row[[joint_col]]]
    }
    
    if (nrow(sub_dt) < 3L || nrow(std_sub) == 0L) {
      return(list(
        fit = NULL,
        diag = data.table::data.table(
          .grp_id = gid, selected_model = "failed_insufficient_data",
          reason = "insufficient_data", total_deaths = total_deaths,
          total_exposure = total_exp, bic = NA_real_, edf = NA_real_,
          converged = FALSE, monotone_ok = FALSE, n_violations = NA_integer_,
          overridden = FALSE
        )
      ))
    }
    
    k_flex <- min(as.integer(dof), max(3L, nrow(sub_dt) - 2L))
    k_low  <- min(5L, k_flex)
    
    # Fitters
    fit_brass <- function(d) {
      tryCatch(
        mgcv::gam(Dx ~ 1 + log_mx_std + offset(log(Ex)),
                  family = mgcv::nb(), data = d, method = "REML"),
        error = function(e) NULL
      )
    }
    fit_flex <- function(d, k_val = k_flex, sel = TRUE) {
      tryCatch(
        mgcv::gam(Dx ~ 1 + log_mx_std + s(x, k = k_val, bs = basis) + offset(log(Ex)),
                  family = mgcv::nb(), data = d, method = "REML", select = sel),
        error = function(e) NULL
      )
    }
    fit_monotonic <- function(d) {
      if (requireNamespace("scam", quietly = TRUE)) {
        m <- tryCatch(
          suppressWarnings(
            scam::scam(Dx ~ 1 + log_mx_std + s(x, k = k_low, bs = "mpi") + offset(log(Ex)),
                       family = stats::quasipoisson(), data = d)
          ),
          error = function(e) NULL
        )
        if (!is.null(m)) return(m)
      }
      fit_flex(d, k_val = k_low, sel = TRUE)
    }
    fit_joint <- function() {
      if (!has_joint) return(NULL)
      curr_val <- g_row[[joint_col]]
      if (curr_val %in% joint_exclude) return(NULL)
      jid <- g_row$.joint_id
      sib_ids <- group_dict[.joint_id == jid & !(get(joint_col) %in% joint_exclude)]$.grp_id
      if (length(sib_ids) < 2L) return(NULL)
      j_dt <- dt[.grp_id %in% sib_ids & .use_fit == TRUE]
      if (nrow(j_dt) < 6L) return(NULL)
      j_dt[, (joint_col) := as.factor(get(joint_col))]
      fml <- stats::as.formula(
        paste0("Dx ~ ", joint_col, " + log_mx_std + s(x, k = ", k_flex,
               ", bs = '", basis, "') + offset(log(Ex))")
      )
      tryCatch(
        mgcv::gam(fml, family = mgcv::nb(), data = j_dt, method = "REML", select = TRUE),
        error = function(e) NULL
      )
    }
    
    # --- PATH 1: MANUAL OVERRIDE ---
    fm_target <- g_row$forced_model
    if (!is.na(fm_target)) {
      res <- switch(
        fm_target,
        "garm_flexible"  = eval_candidate(fit_flex(sub_dt, k_flex, sel = FALSE), std_sub, "garm_flexible"),
        "garm_monotonic" = eval_candidate(fit_monotonic(sub_dt), std_sub, "garm_monotonic", kannisto_tail = TRUE),
        "gam_kannisto"   = eval_candidate(fit_flex(sub_dt, k_flex, sel = TRUE), std_sub, "gam_kannisto", kannisto_tail = TRUE),
        "garm_joint"     = eval_candidate(fit_joint(), std_sub, "garm_joint"),
        "brass"          = eval_candidate(fit_brass(sub_dt), std_sub, "brass"),
        stop("Unknown `forced_model`: ", fm_target)
      )
      return(list(
        fit = res$fit,
        diag = data.table::data.table(
          .grp_id = gid, selected_model = res$model, reason = "manual_override",
          total_deaths = total_deaths, total_exposure = total_exp,
          bic = res$bic, edf = res$edf, converged = res$converged,
          monotone_ok = res$monotone_ok, n_violations = res$n_violations,
          overridden = TRUE
        )
      ))
    }
    
    # --- PATH 2: AUTOMATIC SELECTION + PLAUSIBILITY LADDER ---
    res_brass <- eval_candidate(fit_brass(sub_dt), std_sub, "brass")
    
    if (total_deaths < min_deaths_for_smooth) {
      return(list(
        fit = res_brass$fit,
        diag = data.table::data.table(
          .grp_id = gid, selected_model = "brass", reason = "low_deaths_fallback",
          total_deaths = total_deaths, total_exposure = total_exp,
          bic = res_brass$bic, edf = res_brass$edf, converged = res_brass$converged,
          monotone_ok = res_brass$monotone_ok, n_violations = res_brass$n_violations,
          overridden = FALSE
        )
      ))
    }
    
    # Stage A: Evaluate true Negative Binomial candidates on comparable BIC
    nb_candidates <- list(
      garm_flexible = eval_candidate(fit_flex(sub_dt, k_flex, sel = TRUE), std_sub, "garm_flexible"),
      gam_kannisto  = eval_candidate(fit_flex(sub_dt, k_flex, sel = TRUE), std_sub, "gam_kannisto", kannisto_tail = TRUE),
      brass         = res_brass
    )
    
    if (has_joint && !(g_row[[joint_col]] %in% joint_exclude)) {
      j_fit <- fit_joint()
      if (!is.null(j_fit)) {
        std_j <- data.table::copy(std_sub)
        std_j[, (joint_col) := factor(get(joint_col), levels = levels(j_fit$model[[joint_col]]))]
        nb_candidates$garm_joint <- eval_candidate(j_fit, std_j, "garm_joint")
      }
    }
    
    # Check if any smooth NB model converged and passed 50+ monotonicity
    smooth_nb <- setdiff(names(nb_candidates), "brass")
    valid_smooth <- smooth_nb[vapply(nb_candidates[smooth_nb], function(x) isTRUE(x$ok), logical(1))]
    
    if (length(valid_smooth) > 0L) {
      valid_all <- names(nb_candidates)[vapply(nb_candidates, function(x) isTRUE(x$ok), logical(1))]
      bics <- vapply(valid_all, function(nm) {
        b <- nb_candidates[[nm]]$bic
        if (is.finite(b)) b else Inf
      }, numeric(1))
      best_name <- valid_all[which.min(bics)]
      best <- nb_candidates[[best_name]]
      reason <- if (best_name == "garm_flexible") "bic_and_checks_ok" else "fallback_checks_or_bic"
    } else {
      # Stage B: Unconstrained smooth NB models failed 50+ monotonicity;
      # try shape-restricted monotonic spline (with Kannisto tail) vs. Brass
      res_mono <- eval_candidate(fit_monotonic(sub_dt), std_sub, "garm_monotonic", kannisto_tail = TRUE)
      if (isTRUE(res_mono$ok)) {
        best <- res_mono
        reason <- "monotonic_spline_fallback"
      } else if (isTRUE(res_brass$ok)) {
        best <- res_brass
        reason <- "brass_monotonic_fallback"
      } else {
        best <- res_brass
        reason <- "all_checks_failed_brass_fallback"
      }
    }
    
    list(
      fit = best$fit,
      diag = data.table::data.table(
        .grp_id = gid, selected_model = best$model, reason = reason,
        total_deaths = total_deaths, total_exposure = total_exp,
        bic = best$bic, edf = best$edf, converged = best$converged,
        monotone_ok = best$monotone_ok, n_violations = best$n_violations,
        overridden = FALSE
      )
    )
  }
  
  # Execute across groups (cross-platform parallel or sequential)
  gids <- group_dict$.grp_id
  if (workers > 1L && requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = as.integer(workers))
    results <- future.apply::future_lapply(gids, fit_one_group, future.seed = TRUE)
  } else {
    results <- lapply(gids, fit_one_group)
  }
  
  fits_list <- lapply(results, `[[`, "fit")
  names(fits_list) <- as.character(gids)
  diag_dt <- data.table::rbindlist(lapply(results, `[[`, "diag"), use.names = TRUE)
  diag_dt <- merge(group_dict[, !c("forced_model", ".joint_id"), with = FALSE],
                   diag_dt, by = ".grp_id")
  
  structure(
    list(
      fits = fits_list,
      diagnostics = diag_dt,
      group_dict = group_dict[, !c("forced_model", ".joint_id"), with = FALSE],
      meta = list(
        by = by,
        std_by = std_by,
        pred_ages = pred_ages,
        monotone_age = monotone_age,
        joint_col = joint_col
      )
    ),
    class = "garmort_fits"
  )
}


#' Predict Mortality Rates and Posterior Simulations from Fitted Models
#'
#' @param object A `"garmort_fits"` object returned by `fit_mortality_models()`.
#' @param standard The reference standard `data.table` containing `x` and `log_mx_std`.
#' @param pred_ages Numeric vector of ages to predict (defaults to `object$meta$pred_ages`).
#' @param n_sim Integer number of posterior draws per group via `lpmatrix` and
#'   `MASS::mvrnorm()` (default `0L` for point estimates only).
#' @param seed Optional integer seed for reproducible simulations.
#' @return A list with:
#'   - `estimates`: `data.table` with grouping columns, `x`, `mx_std`, and `mx_s`.
#'   - `simulations`: `data.table` with grouping columns, `x`, and `sim_1..sim_n`
#'     (or `NULL` if `n_sim == 0`).
#' @importFrom data.table as.data.table copy rbindlist setnames fcase
#' @importFrom MASS mvrnorm
#' @importFrom stats coef predict plogis
#' @export
predict_mortality_rates <- function(object,
                                    standard,
                                    pred_ages = NULL,
                                    n_sim = 0L,
                                    seed = NULL) {
  if (!inherits(object, "garmort_fits")) {
    stop("`object` must be a `garmort_fits` object.")
  }
  if (!is.null(seed)) set.seed(seed)
  
  std <- data.table::copy(data.table::as.data.table(standard))
  if (!"log_mx_std" %in% names(std)) std[, log_mx_std := log(mx_std)]
  if (!"mx_std" %in% names(std)) std[, mx_std := exp(log_mx_std)]
  
  ages <- if (is.null(pred_ages)) object$meta$pred_ages else pred_ages
  by_cols <- object$meta$by
  std_by <- object$meta$std_by
  joint_col <- object$meta$joint_col
  gdict <- object$group_dict
  
  est_list <- vector("list", nrow(gdict))
  sim_list <- if (n_sim > 0L) vector("list", nrow(gdict)) else NULL
  
  for (i in seq_len(nrow(gdict))) {
    gid <- gdict$.grp_id[i]
    fit_wrap <- object$fits[[as.character(gid)]]
    if (is.null(fit_wrap)) next
    
    g_row <- gdict[i]
    if (length(std_by) > 0L) {
      nd <- merge(
        data.table::data.table(x = ages),
        std[g_row[, std_by, with = FALSE], on = std_by, nomatch = NULL],
        by = "x"
      )
    } else {
      nd <- std[x %in% ages]
    }
    nd[, Ex := 1]
    
    eng <- fit_wrap$engine
    if (!is.null(joint_col) && joint_col %in% names(g_row)) {
      if (joint_col %in% names(eng$model) && is.factor(eng$model[[joint_col]])) {
        nd[, (joint_col) := factor(g_row[[joint_col]], levels = levels(eng$model[[joint_col]]))]
      } else {
        nd[, (joint_col) := g_row[[joint_col]]]
      }
    }
    
    # Point estimates via lpmatrix (cleanly excludes offset)
    Xp <- stats::predict(eng, newdata = nd, type = "lpmatrix")
    b_raw <- as.numeric(if (!is.null(stats::coef(eng))) stats::coef(eng) else eng$coefficients)
    if (length(b_raw) == 0L && !is.null(eng$coefficients)) {
      b_raw <- as.numeric(eng$coefficients)
    }
    
    is_scam <- inherits(eng, "scam") && !is.null(eng$p.ident)
    b_hat <- b_raw
    if (is_scam) {
      p_id <- eng$p.ident
      b_hat[p_id] <- exp(b_hat[p_id])
    }
    eta_hat <- as.numeric(Xp %*% b_hat)
    mx_hat <- exp(eta_hat)
    
    # Apply Kannisto blending if selected
    hi_idx <- which(nd$x >= 80)
    w_k <- NULL
    mx_k <- NULL
    if (!is.null(fit_wrap$kannisto_mod) && length(hi_idx) > 0L) {
      mx_k <- as.numeric(stats::plogis(stats::predict(fit_wrap$kannisto_mod, newdata = nd[hi_idx])))
      w_k <- data.table::fcase(
        nd$x[hi_idx] < 100, 0,
        nd$x[hi_idx] >= 100 & nd$x[hi_idx] < 105, (nd$x[hi_idx] - 100) / 5,
        default = 1
      )
      mx_hat[hi_idx] <- (1 - w_k) * mx_hat[hi_idx] + w_k * mx_k
    }
    
    est_dt <- data.table::data.table(x = nd$x, mx_std = nd$mx_std, mx_s = mx_hat)
    if (!is.null(by_cols) && length(by_cols) > 0L) {
      for (col in by_cols) est_dt[, (col) := g_row[[col]]]
    }
    est_list[[i]] <- est_dt
    
    # Posterior simulations
    if (n_sim > 0L && length(b_raw) > 0L) {
      if (is_scam && !is.null(eng$Vp.t) && all(dim(eng$Vp.t) == c(length(b_hat), length(b_hat)))) {
        # Simulate directly on delta-method transformed scale for scam
        beta_sims <- MASS::mvrnorm(n_sim, mu = b_hat, Sigma = eng$Vp.t)
        if (is.null(dim(beta_sims))) beta_sims <- matrix(beta_sims, nrow = 1L)
        p_id <- eng$p.ident
        beta_sims[, p_id] <- pmax(beta_sims[, p_id, drop = FALSE], 0)
      } else if (!is.null(eng$Vp) && all(dim(eng$Vp) == c(length(b_raw), length(b_raw)))) {
        beta_sims <- MASS::mvrnorm(n_sim, mu = b_raw, Sigma = eng$Vp)
        if (is.null(dim(beta_sims))) beta_sims <- matrix(beta_sims, nrow = 1L)
        if (is_scam) {
          p_id <- eng$p.ident
          beta_sims[, p_id] <- exp(pmin(beta_sims[, p_id, drop = FALSE], 5))
        }
      } else {
        beta_sims <- NULL
      }
      
      if (!is.null(beta_sims)) {
        eta_sims <- Xp %*% t(beta_sims)
        eta_sims <- pmin(pmax(eta_sims, -30), log(1.5))
        rate_sims <- exp(eta_sims)
        
        if (!is.null(fit_wrap$kannisto_mod) && length(hi_idx) > 0L) {
          rate_sims[hi_idx, ] <- (1 - w_k) * rate_sims[hi_idx, , drop = FALSE] + w_k * mx_k
        }
        
        s_dt <- data.table::as.data.table(rate_sims)
        data.table::setnames(s_dt, paste0("sim_", seq_len(n_sim)))
        s_dt <- cbind(data.table::data.table(x = nd$x), s_dt)
        if (!is.null(by_cols) && length(by_cols) > 0L) {
          for (col in by_cols) s_dt[, (col) := g_row[[col]]]
        }
        sim_list[[i]] <- s_dt
      }
    }
  }
  
  est_out <- data.table::rbindlist(est_list, use.names = TRUE)
  sim_out <- if (n_sim > 0L) data.table::rbindlist(sim_list, use.names = TRUE) else NULL
  
  if (!is.null(by_cols) && length(by_cols) > 0L) {
    data.table::setcolorder(est_out, c(by_cols, "x", "mx_std", "mx_s"))
    if (!is.null(sim_out)) {
      data.table::setcolorder(sim_out, c(by_cols, "x", paste0("sim_", seq_len(n_sim))))
    }
  }
  
  list(estimates = est_out, simulations = sim_out)
}