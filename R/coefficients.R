#' Second-order coefficient of the scaled risk difference
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability, or ability estimate.
#' @param estimator Required: exactly `"MLE"` or `"WLE"`.
#' @return The numeric coefficient Delta_2, based on design averages.
#'   A single item is valid; no response enumeration is needed.
#' @export
#' @examples
#' delta2(1, 0, -1, "MLE")
delta2 <- function(a, b, theta, estimator) {
  validate_design(a, b, theta)
  validate_estimator(estimator)
  k <- coefficient_cumulants(a, b, theta)
  unname(k[["k3"]] * delta2_factor(k, estimator))
}

coefficient_cumulants <- function(a, b, theta) {
  k <- design_cumulants(a, b, theta)
  if (any(!is.finite(k)) || k[["k2"]] <= 0) {
    stop("a, b and theta yield non-finite cumulants or zero information.",
         call. = FALSE)
  }
  k
}

delta2_factor <- function(k, estimator) {
  k2 <- k[["k2"]]
  k3 <- k[["k3"]]
  k4 <- k[["k4"]]
  k5 <- k[["k5"]]
  if (estimator == "MLE") {
    -(35 * k3^3 - 30 * k2 * k3 * k4 + 4 * k2^2 * k5) / (8 * k2^7)
  } else {
    k3 * (5 * k2 * k4 - 8 * k3^2) / (8 * k2^7)
  }
}

#' Third-order coefficient of the scaled risk difference
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability, or ability estimate.
#' @param estimator Required: exactly `"MLE"` or `"WLE"`.
#' @return The numeric coefficient Delta_3 for the MLE. For the WLE,
#'   `NA_real_` and a message stating that the coefficient is unavailable.
#'   A single item is valid; no response enumeration is needed.
#' @export
#' @examples
#' delta3(1, 0, 0, "MLE")
delta3 <- function(a, b, theta, estimator) {
  validate_design(a, b, theta)
  validate_estimator(estimator)
  if (estimator == "WLE") {
    message("No third-order coefficient is available for the WLE.")
    return(NA_real_)
  }
  k <- coefficient_cumulants(a, b, theta)
  delta3_mle(k)
}

delta3_mle <- function(k) {
  k2 <- k[["k2"]]
  k3 <- k[["k3"]]
  k4 <- k[["k4"]]
  k5 <- k[["k5"]]
  k6 <- k[["k6"]]
  k7 <- k[["k7"]]
  r3 <- 288 * k2^4 * k3 * k7 + 288 * k2^4 * k4 * k6 +
    384 * k2^4 * k5^2 - 5152 * k2^3 * k3^2 * k6 -
    16928 * k2^3 * k3 * k4 * k5 - 2120 * k2^3 * k4^3 +
    49552 * k2^2 * k3^3 * k5 + 93388 * k2^2 * k3^2 * k4^2 -
    279274 * k2 * k3^4 * k4 + 169263 * k3^6
  -r3 / (768 * k2^10)
}
