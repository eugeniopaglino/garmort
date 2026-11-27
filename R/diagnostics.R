#' Plot Observed vs. Modeled vs. Standard Log-Mortality Rates
#'
#' @param estimates `data.table` from `predict_mortality_rates()$estimates`.
#' @param observed Optional `data.table` of observed rates (with `x` and `mx` or `Dx`/`Ex`).
#' @param facet_formula Optional formula passed to `ggplot2::facet_grid()` or `facet_wrap()`.
#' @return A `ggplot` object.
#' @export
plot_mortality_fit <- function(estimates, observed = NULL, facet_formula = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for `plot_mortality_fit()`.")
  }
  est <- data.table::as.data.table(estimates)
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = est, ggplot2::aes(x = x, y = log(mx_s), color = "Modeled"), linewidth = 1) +
    ggplot2::geom_line(data = est, ggplot2::aes(x = x, y = log(mx_std), color = "Standard"),
                       linewidth = 0.8, linetype = "dashed")
  
  if (!is.null(observed)) {
    obs <- data.table::copy(data.table::as.data.table(observed))
    if (!"mx" %in% names(obs) && all(c("Dx", "Ex") %in% names(obs))) {
      obs[, mx := Dx / Ex]
    }
    p <- p + ggplot2::geom_point(
      data = obs[mx > 0],
      ggplot2::aes(x = x, y = log(mx), color = "Observed"),
      alpha = 0.7, size = 1.5
    )
  }
  
  if (!is.null(facet_formula)) {
    p <- p + ggplot2::facet_wrap(facet_formula)
  }
  
  p +
    ggplot2::labs(x = "Age", y = "Log Mortality Rate", color = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}