#' Compute Single-Year Life Tables (Grouped or Single Series)
#'
#' @param data A `data.table` (or data.frame) with age and mortality rate columns.
#' @param age_col Name of the exact integer age column (default `"x"`).
#' @param mx_col Name of the mortality rate column (default `"mx_s"`).
#' @param by Optional character vector of grouping columns.
#' @param radix Life table radix (default `1e5`).
#' @param ax_method Character: `"keyfitz"` (Keyfitz-Preston infant rule `0.07 + 1.7 * mx[1]`,
#'   `0.5` for intermediate ages, `1/mx[N]` for open interval) or `"exponential"`
#'   (uses `solve_ax_exponential()` from local log-mortality slopes).
#' @return A `data.table` with `by` columns (if any), `age_col`, and
#'   `mx, ax, qx, px, lx, dx, Lx, Tx, ex`.
#' @importFrom data.table as.data.table copy setnames
#' @export
lt_single <- function(data,
                      age_col = "x",
                      mx_col = "mx_s",
                      by = NULL,
                      radix = 1e5,
                      ax_method = c("keyfitz", "exponential")) {
  ax_method <- match.arg(ax_method)
  dt <- data.table::copy(data.table::as.data.table(data))
  if (!age_col %in% names(dt)) stop("Column `", age_col, "` not found.")
  if (!mx_col %in% names(dt)) stop("Column `", mx_col, "` not found.")
  
  calc_one <- function(age, mx) {
    ord <- order(age)
    age <- as.numeric(age[ord])
    mx  <- as.numeric(mx[ord])
    N   <- length(mx)
    
    if (ax_method == "keyfitz") {
      ax <- rep(0.5, N)
      if (age[1] == 0) ax[1] <- pmin(0.07 + 1.7 * mx[1], 0.5)
      ax[N] <- 1 / mx[N]
    } else {
      log_mx <- log(pmax(mx, 1e-12))
      beta <- c(diff(log_mx), log_mx[N] - log_mx[N - 1L])
      ax <- solve_ax_exponential(n = rep(1, N), beta = beta, r = 0)
      if (age[1] == 0) ax[1] <- pmin(0.07 + 1.7 * mx[1], 0.5)
      ax[N] <- 1 / mx[N]
    }
    
    qx <- mx / (1 + (1 - ax) * mx)
    qx[N] <- 1
    qx <- pmin(pmax(qx, 0), 1)
    px <- 1 - qx
    lx <- c(radix, radix * cumprod(px[-N]))
    dx <- -diff(c(lx, 0))
    Lx <- c(lx[-1] + ax[-N] * dx[-N], lx[N] * ax[N])
    Tx <- rev(cumsum(rev(Lx)))
    ex <- Tx / lx
    
    list(.age_out = age, mx = mx, ax = ax, qx = qx, px = px,
         lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex)
  }
  
  out <- if (is.null(by) || length(by) == 0L) {
    dt[, calc_one(get(age_col), get(mx_col))]
  } else {
    dt[, calc_one(get(age_col), get(mx_col)), by = by]
  }
  data.table::setnames(out, ".age_out", age_col)
  out
}


#' Compute Abridged Life Tables (Variable Start/End Ages, Grouped or Single Series)
#'
#' @param data A `data.table` (or data.frame) with interval starting ages and rates.
#' @param age_col Name of the interval starting age column (default `"x"`).
#' @param mx_col Name of the mortality rate column (default `"nMx"`).
#' @param by Optional character vector of grouping columns.
#' @param radix Life table radix (default `1e5`).
#' @param ax_method Character: `"standard"` (Keyfitz-Preston for age 0, 0.4 for 1-4,
#'   0.3 for 0-4, 0.5 otherwise, `1/nMx[N]` for open interval) or `"exponential"`.
#' @return A `data.table` with `by` columns (if any), `age_col`, and
#'   `n, nMx, nax, nqx, npx, lx, ndx, nLx, Tx, ex`.
#' @importFrom data.table as.data.table copy setnames
#' @export
lt_abridged <- function(data,
                        age_col = "x",
                        mx_col = "nMx",
                        by = NULL,
                        radix = 1e5,
                        ax_method = c("standard", "exponential")) {
  ax_method <- match.arg(ax_method)
  dt <- data.table::copy(data.table::as.data.table(data))
  if (!age_col %in% names(dt)) stop("Column `", age_col, "` not found.")
  if (!mx_col %in% names(dt)) stop("Column `", mx_col, "` not found.")
  
  calc_one <- function(age, nMx) {
    ord <- order(age)
    age <- as.numeric(age[ord])
    nMx <- as.numeric(nMx[ord])
    N   <- length(nMx)
    w   <- c(diff(age), 5)
    
    if (ax_method == "standard") {
      ax_frac <- rep(0.5, N)
      if (age[1] == 0 && w[1] == 1) {
        ax_frac[1] <- pmin(0.07 + 1.7 * nMx[1], 0.5)
        if (N >= 2L && age[2] == 1 && w[2] == 4) ax_frac[2] <- 0.4
      } else if (age[1] == 0 && w[1] == 5) {
        ax_frac[1] <- 0.3
      }
      nax <- w * ax_frac
      nax[N] <- 1 / nMx[N]
    } else {
      log_mx <- log(pmax(nMx, 1e-12))
      beta <- c(diff(log_mx) / w[-N], (log_mx[N] - log_mx[N - 1L]) / w[N - 1L])
      nax <- solve_ax_exponential(n = w, beta = beta, r = 0)
      if (age[1] == 0 && w[1] == 1) nax[1] <- pmin(0.07 + 1.7 * nMx[1], 0.5)
      nax[N] <- 1 / nMx[N]
    }
    
    nqx <- (w * nMx) / (1 + (w - nax) * nMx)
    nqx[N] <- 1
    nqx <- pmin(pmax(nqx, 0), 1)
    npx <- 1 - nqx
    lx  <- c(radix, radix * cumprod(npx[-N]))
    ndx <- -diff(c(lx, 0))
    nLx <- c(w[-N] * lx[-1] + nax[-N] * ndx[-N], lx[N] * nax[N])
    Tx  <- rev(cumsum(rev(nLx)))
    ex  <- Tx / lx
    
    list(.age_out = age, n = w, nMx = nMx, nax = nax, nqx = nqx, npx = npx,
         lx = lx, ndx = ndx, nLx = nLx, Tx = Tx, ex = ex)
  }
  
  out <- if (is.null(by) || length(by) == 0L) {
    dt[, calc_one(get(age_col), get(mx_col))]
  } else {
    dt[, calc_one(get(age_col), get(mx_col)), by = by]
  }
  data.table::setnames(out, ".age_out", age_col)
  out
}


#' Aggregate Posterior Rate Simulations into Single-Year & Abridged Life Table Summaries
#'
#' @param simulations A `data.table` from `predict_mortality_rates(n_sim > 0)$simulations`.
#' @param age_col Name of the age column in `simulations` (default `"x"`). Preserved in output.
#' @param by Optional character vector of grouping columns.
#' @param max_abridged_age Starting age of the open-ended interval for abridged
#'   summaries (default `85`).
#' @param radix Life table radix (default `1`).
#' @param probs Numeric length-2 vector of credible interval quantiles (default `c(0.025, 0.975)`).
#' @return A list with `single_year_summary` and `abridged_summary` `data.table`s.
#' @importFrom data.table as.data.table copy data.table rbindlist setnames setorderv setcolorder fifelse
#' @importFrom stats quantile
#' @export
sim_life_table <- function(simulations,
                           age_col = "x",
                           by = NULL,
                           max_abridged_age = 85,
                           radix = 1,
                           probs = c(0.025, 0.975)) {
  sim_dt <- data.table::copy(data.table::as.data.table(simulations))
  if (!age_col %in% names(sim_dt)) stop("Column `", age_col, "` not found in `simulations`.")
  
  sim_cols <- grep("^sim_", names(sim_dt), value = TRUE)
  if (length(sim_cols) == 0L) stop("No `sim_*` columns found in `simulations`.")
  
  if (is.null(by) || length(by) == 0L) {
    sim_dt[, .grp_id := 1L]
    gdict <- data.table::data.table(.grp_id = 1L)
  } else {
    sim_dt[, .grp_id := .GRP, by = by]
    gdict <- unique(sim_dt[, c(".grp_id", by), with = FALSE])
  }
  
  summarize_one <- function(sub_dt) {
    data.table::setorderv(sub_dt, age_col)
    exact_ages <- as.integer(floor(sub_dt[[age_col]]))
    n_ages <- length(exact_ages)
    n_sims <- length(sim_cols)
    
    min_age <- min(exact_ages)
    abridged_starts <- if (min_age == 0L) {
      c(0L, 1L, seq(5L, as.integer(max_abridged_age), by = 5L))
    } else if (min_age == 1L) {
      c(1L, seq(5L, as.integer(max_abridged_age), by = 5L))
    } else {
      seq(min_age, as.integer(max_abridged_age), by = 5L)
    }
    
    age_groups <- cut(exact_ages, breaks = c(abridged_starts, Inf),
                      right = FALSE, labels = abridged_starts)
    
    s_metrics <- c("mx", "lx", "ex")
    a_metrics <- c("nMx", "nax", "nqx", "npx", "ndx", "lx", "nLx", "Tx", "ex")
    
    s_mats <- stats::setNames(lapply(s_metrics, function(m) matrix(NA_real_, n_ages, n_sims)), s_metrics)
    a_mats <- stats::setNames(lapply(a_metrics, function(m) matrix(NA_real_, length(abridged_starts), n_sims)), a_metrics)
    
    for (j in seq_len(n_sims)) {
      mx_vec <- sub_dt[[sim_cols[j]]]
      slt <- lt_single(data.table::data.table(x = exact_ages, mx_s = mx_vec),
                       age_col = "x", mx_col = "mx_s", radix = radix)
      
      for (m in s_metrics) s_mats[[m]][, j] <- slt[[m]]
      
      slt[, age_group := age_groups]
      agg <- slt[, .(nLx = sum(Lx), ndx = sum(dx), lx = lx[1], Tx = Tx[1], ex = ex[1], n = .N), by = age_group]
      agg[, nMx := ndx / nLx]
      agg[, nqx := ndx / lx]
      agg[.N, nqx := 1]
      agg[, npx := 1 - nqx]
      agg[, nax := data.table::fifelse(ndx > 0, (nLx - n * (lx - ndx)) / ndx, n / 2)]
      agg[.N, nax := ex]
      
      for (m in a_metrics) a_mats[[m]][, j] <- agg[[m]]
    }
    
    collapse_mats <- function(mat_list, ages_vec) {
      parts <- lapply(names(mat_list), function(m) {
        mat <- mat_list[[m]]
        out <- data.table::data.table(
          mean = rowMeans(mat, na.rm = TRUE),
          l95  = apply(mat, 1, stats::quantile, probs = probs[1], na.rm = TRUE),
          u95  = apply(mat, 1, stats::quantile, probs = probs[2], na.rm = TRUE)
        )
        data.table::setnames(out, paste0(m, c("_mean", "_l95", "_u95")))
        out
      })
      age_dt <- data.table::data.table(ages_vec)
      data.table::setnames(age_dt, age_col)
      cbind(age_dt, do.call(cbind, parts))
    }
    
    list(
      single = collapse_mats(s_mats, exact_ages),
      abridged = collapse_mats(a_mats, abridged_starts)
    )
  }
  
  s_list <- vector("list", nrow(gdict))
  a_list <- vector("list", nrow(gdict))
  for (i in seq_len(nrow(gdict))) {
    gid <- gdict$.grp_id[i]
    res <- summarize_one(sim_dt[.grp_id == gid])
    if (!is.null(by) && length(by) > 0L) {
      for (col in by) {
        res$single[, (col) := gdict[[col]][i]]
        res$abridged[, (col) := gdict[[col]][i]]
      }
    }
    s_list[[i]] <- res$single
    a_list[[i]] <- res$abridged
  }
  
  s_out <- data.table::rbindlist(s_list, use.names = TRUE)
  a_out <- data.table::rbindlist(a_list, use.names = TRUE)
  
  if (!is.null(by) && length(by) > 0L) {
    data.table::setcolorder(s_out, c(by, age_col))
    data.table::setcolorder(a_out, c(by, age_col))
  }
  
  list(
    single_year_summary = s_out,
    abridged_summary    = a_out
  )
}