#' Second-order coefficient of the scaled risk difference
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability, or ability estimate.
#' @param estimator Required: exactly `"MLE"` or `"WLE"`.
#' @return The numeric coefficient Delta_2, based on design averages.
#'   A single item is valid; no response enumeration is needed.
#' @details
#' Cumulant ratios and log information avoid underflow or overflow in powers.
#' A result outside the range of double precision stops with an error containing
#' `representable`. Exact zero coefficients remain zero.
#' @export
#' @examples
#' delta2(1, 0, -1, "MLE")
delta2 <- function(a, b, theta, estimator) {
  validate_design(a, b, theta)
  validate_estimator(estimator)
  design <- design_ratios(a, b, theta)
  polynomial <- delta2_polynomial(design$ratios, estimator)
  exp_value(polynomial[["log_abs"]] - log(8) + 4 * design$log_scale -
              3 * design$log_k2, polynomial[["sign"]], "Delta_2")
}

ratio_factors <- function(r, estimator) {
  powers <- if (estimator == "MLE") {
    rbind(c(3, 0, 0, 0, 0), c(1, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
  } else {
    rbind(c(1, 1, 0, 0, 0), c(3, 0, 0, 0, 0))
  }
  coefficients <- if (estimator == "MLE") c(-35, 30, -4) else c(5, -8)
  second <- polynomial_log(r, powers, coefficients)
  # Signs are positive rescalings of the factors away from their zeros.
  c(sign(r[["k3"]]), second[["sign"]])
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
#' @details
#' Cumulant ratios and log information avoid underflow or overflow in powers.
#' A result outside the range of double precision stops with an error containing
#' `representable`. Exact zero coefficients remain zero.
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
  design <- design_ratios(a, b, theta)
  polynomial <- delta3_polynomial(design$ratios)
  exp_value(polynomial[["log_abs"]] - log(768) + 6 * design$log_scale -
              4 * design$log_k2, -polynomial[["sign"]], "Delta_3")
}

polynomial_log <- function(r, powers, coefficients) {
  r <- r[paste0("k", 3:7)]
  logs <- log(abs(r))
  logs[r == 0] <- 0
  term_logs <- as.vector(powers %*% logs) + log(abs(coefficients))
  signs <- sign(coefficients) * apply(powers, 1L, function(p) prod(sign(r)^p))
  value <- signed_log_mean(term_logs, signs)
  value[["log_abs"]] <- value[["log_abs"]] + log(length(coefficients))
  value
}

delta2_polynomial <- function(r, estimator) {
  powers <- rbind(c(4, 0, 0, 0, 0), c(2, 1, 0, 0, 0))
  if (estimator == "WLE") return(polynomial_log(r, powers, c(-8, 5)))
  powers <- rbind(powers, c(1, 0, 1, 0, 0))
  polynomial_log(r, powers, c(-35, 30, -4))
}

delta3_polynomial <- function(r) {
  powers <- rbind(c(1, 0, 0, 0, 1), c(0, 1, 0, 1, 0), c(0, 0, 2, 0, 0),
                  c(2, 0, 0, 1, 0), c(1, 1, 1, 0, 0), c(0, 3, 0, 0, 0),
                  c(3, 0, 1, 0, 0), c(2, 2, 0, 0, 0), c(4, 1, 0, 0, 0),
                  c(6, 0, 0, 0, 0))
  polynomial_log(r, powers, c(288, 288, 384, -5152, -16928, -2120,
                             49552, 93388, -279274, 169263))
}
