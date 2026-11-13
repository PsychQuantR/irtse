#' Information-based standard error
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability, or ability estimate.
#' @return The actual standard error `1 / sqrt(sum(a^2 * P * (1 - P)))`,
#'   without multiplication by `sqrt(J)`.
#' @export
#' @examples
#' se_info(rep(1, 6), rep(0, 6), 0)
se_info <- function(a, b, theta) {
  validate_design(a, b, theta)
  eta <- a * (theta - b)
  1 / sqrt(sum(a^2 * plogis(eta) * plogis(-eta)))
}

#' Exact parametric-bootstrap standard error
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite generating ability, or plug-in ability estimate.
#' @param estimator Required: exactly `"MLE"` or `"WLE"`.
#' @return The actual standard deviation of the ability estimator at `theta`,
#'   without multiplication by `sqrt(J)`.
#' @details
#' The MLE law is conditional on responses that are neither all correct nor
#' all wrong and requires at least two items. The WLE law includes all responses.
#' Common discriminations use the Lord--Wingersky recursion; commensurable
#' discriminations use lattice convolution; other designs enumerate patterns.
#' Computation stops if the atom count exceeds 100000 or lattice width exceeds
#' 1e6. The package stops instead of approximating. A missing or non-unique
#' estimator root also causes an error naming the sufficient statistic T.
#' @export
#' @examples
#' se_boot(rep(1, 4), rep(0, 4), 0, "MLE")
se_boot <- function(a, b, theta, estimator) {
  validate_design(a, b, theta)
  validate_law_estimator(a, estimator)
  plan <- atom_plan(a, estimator, 100000)
  estimates <- estimate_atoms(plan$T, a, b, estimator)
  law_sd(plan, estimates, a, b, theta)
}
