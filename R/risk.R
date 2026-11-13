#' Exact finite-length scaled risks
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite true generating ability.
#' @param estimator Required: exactly `"MLE"` or `"WLE"`.
#' @param margin Optional single finite non-negative number, with no default.
#'   Only when supplied, also evaluate the rule used by [se_select()].
#' @return A list, in order: `J`, the number of items; `h_true`, the true
#'   bootstrap standard deviation multiplied by `sqrt(J)`; `smse_info` and
#'   `smse_boot`, the mean squared errors of the two SE reports after multiplying
#'   each SE by `sqrt(J)`; and `diff`, their difference. `diff` is not multiplied
#'   by J^2. When `margin` is supplied, append `smse_select`, the scaled risk
#'   of the selection rule, and `p_select`, a numeric vector named `info`,
#'   `boot`, `none` giving the probabilities of those choices. If `none` has
#'   positive probability, `smse_select` is `NA_real_`, including when that
#'   probability underflows to zero in `p_select`. Every atom is checked.
#' @details
#' Both outer expectations and inner bootstrap variances are conditional on
#' non-extreme responses for the MLE; the WLE includes every response. The MLE
#' requires at least two items. Estimates and SE reports are computed once
#' per design and reused. The difference is accumulated in factored form to
#' reduce cancellation. Unlike [se_info()], [se_boot()] and [se_select()],
#' this function uses scaled reports, multiplied by `sqrt(J)`.
#'
#' Exact enumeration stops above 2000 atoms or lattice width 1e6. The package
#' stops instead of approximating. Missing or non-unique roots cause an error
#' naming T. The optional rule is provisional; see [se_select()] for its scope.
#' WLE roots are counted on a grid and at function extrema located by derivative
#' sign changes and searches at sampled local extrema of the derivative.
#' Features narrower than the grid spacing can still be missed; this count
#' does not certify global uniqueness.
#' @template numerical-accuracy
#' @export
#' @examples
#' risk_diff(rep(1, 4), rep(0, 4), 0, "MLE")
risk_diff <- function(a, b, theta, estimator, margin) {
  validate_design(a, b, theta)
  validate_law_estimator(a, estimator)
  if (!missing(margin)) validate_margin(margin)
  eta <- checked_logits(a * (theta - b))
  b <- reference_coordinates(b)$b
  plan <- atom_plan(a, estimator, 2000, "risk_diff()")
  estimates <- estimate_atoms(plan$T, a, b, estimator)
  reports <- scaled_reports(plan, estimates, a, b)
  logs <- atom_log_probabilities(plan, eta)
  weights <- exp(logs)
  J <- length(a)
  h <- sqrt(J) * weighted_sd(logs, estimates)
  A <- reports$info
  B <- reports$boot
  result <- list(J = J, h_true = h,
                 smse_info = sum(weights * (A - h)^2),
                 smse_boot = sum(weights * (B - h)^2),
                 diff = sum(weights * (A - B) * (A + B - 2 * h)))
  finite_value(unlist(result), "scaled risk")
  if (!missing(margin)) {
    selection <- selection_risk(a, b, estimates, estimator, margin,
                                weights, A, B, h)
    result <- c(result, selection)
  }
  result
}

scaled_reports <- function(plan, estimates, a, b) {
  values <- unique(estimates)
  info <- vapply(values, function(t) sqrt(length(a)) * information_se(a, a * (t - b)),
                 numeric(1))
  boot <- vapply(values, function(t) {
    sqrt(length(a)) * law_sd(plan, estimates, a * (t - b))
  }, numeric(1))
  indices <- match(estimates, values)
  list(info = info[indices], boot = boot[indices])
}

selection_risk <- function(a, b, estimates, estimator, margin, weights, A, B, h) {
  choices <- vapply(estimates, function(t) {
    selection_at(a, a * (t - b), estimator, margin)
  }, character(1))
  probabilities <- vapply(c("info", "boot", "none"), function(choice) {
    sum(weights[choices == choice])
  }, numeric(1))
  report <- ifelse(choices == "info", A, B)
  risk <- if (any(choices == "none")) NA_real_ else {
    finite_value(sum(weights * (report - h)^2), "selection risk")
  }
  list(smse_select = risk, p_select = probabilities)
}
