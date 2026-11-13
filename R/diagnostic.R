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
#' @export
#' @examples
#' se_diagnostic(rep(1, 6), rep(0, 6), -0.5, "MLE")
se_diagnostic <- function(a, b, theta, estimator) {
  validate_design(a, b, theta)
  validate_estimator(estimator)
  d2 <- delta2(a, b, theta, estimator)
  d3 <- if (estimator == "MLE") delta3(a, b, theta, estimator) else NA_real_
  info <- se_info(a, b, theta)
  nearest <- nearest_coefficient_zero(a, b, theta, estimator, info)
  list(delta2 = d2, delta3 = d3, sign_delta2 = as.double(sign(d2)),
       se_info = info, nearest_zero = nearest,
       zero_distance_se = abs(nearest - theta) / info)
}

coefficient_factors <- function(a, b, theta, estimator) {
  eta <- a * (theta - b)
  p <- plogis(eta)
  u <- p * plogis(-eta)
  log_weights <- 2 * log(a) - abs(eta) - 2 * log1p(exp(-abs(eta)))
  if (any(!is.finite(log_weights))) {
    stop("a, b and theta yield non-finite diagnostic factors.", call. = FALSE)
  }
  weights <- exp(log_weights - max(log_weights))
  weights <- weights / sum(weights)
  r3 <- sum(weights * a * (1 - 2 * p))
  r4 <- sum(weights * a^2 * (1 - 6 * u))
  r5 <- sum(weights * a^3 * (1 - 2 * p) * (1 - 12 * u))
  # These are k3 / k2 and 8 * k2^4 * F. Both scale factors are
  # strictly positive, preserving every zero and sign without tiny denominators.
  second <- if (estimator == "MLE") {
    -(35 * r3^3 - 30 * r3 * r4 + 4 * r5)
  } else {
    r3 * (5 * r4 - 8 * r3^2)
  }
  c(r3, second)
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

nearest_coefficient_zero <- function(a, b, theta, estimator, info) {
  grid <- seq(theta - 4 * info, theta + 4 * info, length.out = 4001L)
  factors <- vapply(grid, function(t) {
    coefficient_factors(a, b, t, estimator)
  }, numeric(2))
  candidates <- numeric()
  for (i in seq_len(2L)) {
    f <- function(t) coefficient_factors(a, b, t, estimator)[[i]]
    candidates <- c(candidates, factor_zeros(grid, factors[i, ], f))
  }
  if (!length(candidates)) return(NA_real_)
  candidates[[which.min(abs(candidates - theta))]]
}
