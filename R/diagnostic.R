#' Risk-coefficient diagnostic without choosing a standard error
#'
#' Reports coefficients and the nearest located zero of Delta_2. It does not
#' choose a standard error.
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability estimate.
#' @param estimator Required: exactly `"MLE"` or `"WLE"`.
#' @return A list, in order: `delta2`; `delta3` (`NA_real_` for the WLE,
#'   without a message); `sign_delta2` (a double -1, 0 or 1); `se_info`
#'   (the actual information SE); `nearest_zero`; and `zero_distance_se`
#'   (the distance from `theta` to that zero in information SE units).
#'   The last two elements are `NA_real_` when no zero is located.
#' @details
#' The zero search is limited to the window of 4 information SEs on either side
#' of `theta`. Each of the two factors of Delta_2 is evaluated separately on
#' 4001 equally spaced points. Exact grid zeros and sign changes are candidates;
#' sign changes are refined with absolute root tolerance 1e-12. A factor that
#' is rescaled by a strictly positive quantity has the same zeros and signs;
#' this equivalent rescaling avoids numerical underflow in the search.
#' A factor that
#' touches zero without changing sign can be missed, as can two zeros of the
#' same factor closer together than the grid spacing. This is a local search,
#' not a certificate that no other zeros exist.
#'
#' No estimator-specific response law is used, so no enumeration limit applies
#' and a single item is valid.
#' @template numerical-accuracy
#' @export
#' @examples
#' se_diagnostic(rep(1, 6), rep(0, 6), -0.5, "MLE")
se_diagnostic <- function(a, b, theta, estimator) {
  validate_design(a, b, theta)
  validate_estimator(estimator)
  d2 <- delta2(a, b, theta, estimator)
  d3 <- if (estimator == "MLE") delta3(a, b, theta, estimator) else NA_real_
  info <- se_info(a, b, theta)
  nearest <- nearest_coefficient_zero(a, a * (theta - b), estimator, info)
  zero <- if (is.na(nearest)) NA_real_ else finite_value(theta + nearest, "nearest zero")
  list(delta2 = d2, delta3 = d3, sign_delta2 = as.double(sign(d2)),
       se_info = info, nearest_zero = zero,
       zero_distance_se = abs(nearest) / info)
}

coefficient_factors <- function(a, b, theta, estimator) {
  logit_factors(a, a * (theta - b), estimator)
}

logit_factors <- function(a, eta, estimator) {
  ratio_factors(design_ratios(a, eta, zero_if_unresolved = TRUE)$ratios,
                estimator, zero_if_unresolved = TRUE)
}

factor_zeros <- function(grid, values, f) {
  hits <- grid[values == 0]
  changes <- which(sign(values[-length(values)]) * sign(values[-1L]) < 0)
  refined <- vapply(changes, function(i) {
    stats::uniroot(f, c(grid[[i]], grid[[i + 1L]]),
                   f.lower = values[[i]], f.upper = values[[i + 1L]],
                   tol = 1e-12, check.conv = TRUE)$root
  }, numeric(1))
  c(hits, refined)
}

nearest_coefficient_zero <- function(a, eta, estimator, info) {
  window <- finite_value(c(-4 * info, 4 * info), "diagnostic window")
  grid <- seq(window[[1L]], window[[2L]], length.out = 4001L)
  factors <- vapply(grid, function(d) {
    logit_factors(a, eta + a * d, estimator)
  }, numeric(2))
  candidates <- numeric()
  for (i in seq_len(2L)) {
    f <- function(d) logit_factors(a, eta + a * d, estimator)[[i]]
    candidates <- c(candidates, factor_zeros(grid, factors[i, ], f))
  }
  if (!length(candidates)) return(NA_real_)
  candidates[[which.min(abs(candidates))]]
}
