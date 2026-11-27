#' Compute Interval Widths and Representative Age Midpoints
#'
#' Handles arbitrary starting and ending ages for either abridged or single-year
#' mortality schedules.
#'
#' @param age Numeric vector of interval starting ages (sorted ascending).
#' @param open_width Numeric width assumed for the final open-ended interval
#'   when computing a representative plotting/matching midpoint (default 5).
#' @param infant_mid Numeric representative age for the [0, 1) interval (default 0.2).
#' @param child04_mid Numeric representative age for a combined [0, 5) interval (default 1.5).
#' @return A `data.table` with columns `age`, `n` (width), `x_mid`, and `is_open`.
#' @importFrom data.table data.table fifelse
#' @export
compute_age_intervals <- function(age,
                                  open_width = 5,
                                  infant_mid = 0.2,
                                  child04_mid = 1.5) {
  age <- sort(unique(as.numeric(age)))
  k <- length(age)
  if (k < 2L) stop("At least two age groups are required.")
  
  n <- c(diff(age), open_width)
  is_open <- c(rep(FALSE, k - 1L), TRUE)
  
  x_mid <- age + n / 2
  # Special handling for infant [0, 1) and combined [0, 5)
  x_mid[age == 0 & n == 1] <- infant_mid
  x_mid[age == 0 & n == 5] <- child04_mid
  
  data.table::data.table(
    age = age,
    n = n,
    x_mid = x_mid,
    is_open = is_open
  )
}

#' Analytical Fraction of Interval Lived (a_x) Under Exponential Hazard & Exposure
#'
#' Computes the implied average age-offset `a_x` within an interval `[x, x + n)`
#' assuming log-linear hazard `mu(x + t) = mu0 * exp(beta * t)` and exponential
#' exposure density `N(x + t) = N0 * exp(r * t)`.
#'
#' @param n Numeric vector of interval widths.
#' @param beta Numeric vector of log-mortality slopes within each interval.
#' @param r Numeric vector of log-exposure slopes within each interval (default 0).
#' @return Numeric vector of `a_x` values (in years, between `0` and `n`).
#' @export
solve_ax_exponential <- function(n, beta, r = 0) {
  s <- beta + r
  ax <- ifelse(
    abs(s * n) < 1e-6,
    n / 2 + (s * n^2) / 12,
    (n * exp(s * n)) / (exp(s * n) - 1) - (1 / s)
  )
  pmin(pmax(ax, 0), n)
}

#' Check Non-Decreasing Mortality Past a Threshold Age
#'
#' @param x Numeric vector of ages (will be sorted).
#' @param mx Numeric vector of mortality rates (or log-rates).
#' @param min_age Numeric threshold age above which mortality must not decline (default 50).
#' @param tol Numeric tolerance for tiny numerical noise (default 0).
#' @return A list with `ok` (logical), `n_violations` (integer), and `max_drop` (numeric).
#' @export
check_monotonicity <- function(x, mx, min_age = 50, tol = 0) {
  ord <- order(x)
  x <- x[ord]
  mx <- mx[ord]
  
  keep <- which(x >= min_age & is.finite(mx))
  if (length(keep) < 2L) {
    return(list(ok = TRUE, n_violations = 0L, max_drop = 0))
  }
  
  diffs <- diff(mx[keep])
  violations <- which(diffs < -tol)
  list(
    ok = length(violations) == 0L,
    n_violations = length(violations),
    max_drop = if (length(violations) > 0L) min(diffs[violations]) else 0
  )
}