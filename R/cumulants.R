#' Design averages of logistic cumulants
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability, or ability estimate.
#' @return A numeric vector named `k2`, ..., `k7`, with element r equal to
#'   the mean of `a^r` times the rth Bernoulli cumulant. A single item is valid.
#' @export
#' @examples
#' design_cumulants(1, 0, 0)
design_cumulants <- function(a, b, theta) {
  validate_design(a, b, theta)
  eta <- a * (theta - b)
  p <- plogis(eta)
  u <- p * plogis(-eta)
  s <- 1 - 2 * p
  c(k2 = mean(a^2 * u),
    k3 = mean(a^3 * u * s),
    k4 = mean(a^4 * u * (1 - 6 * u)),
    k5 = mean(a^5 * u * s * (1 - 12 * u)),
    k6 = mean(a^6 * u * (1 - 30 * u + 120 * u^2)),
    k7 = mean(a^7 * u * s * (1 - 60 * u + 360 * u^2)))
}
