#' Provisional choice between two standard errors
#'
#' This rule is not a recommendation. Its form is provisional.
#' It carries no guarantee at a fixed test length
#' or near a zero of Delta_2, and `margin` has no theoretical default.
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability estimate at which to evaluate the rule.
#' @param estimator Required: exactly `"MLE"` or `"WLE"`.
#' @param margin Required single finite non-negative number. The threshold
#'   applies to the absolute value of Delta_2; it has no theoretical default.
#' @return A list in the following order:
#' \describe{
#'   \item{choice}{`"info"`, `"boot"`, or `"none"` when no choice is made.}
#'   \item{se}{The chosen actual standard error, or `NA_real_` for `"none"`.}
#'   \item{se_info}{The actual information standard error.}
#'   \item{se_boot}{The actual exact bootstrap standard error.}
#'   \item{delta2}{The second-order scaled risk coefficient.}
#'   \item{delta3}{The third-order coefficient for the MLE, or `NA_real_`
#'     for the WLE, without emitting a message.}
#'   \item{within_margin}{Whether `abs(delta2) <= margin`.}
#'   \item{reason}{One sentence identifying the deciding coefficient and sign,
#'     or explaining why no choice is made.}
#' }
#' The three SE elements are actual standard errors, without `sqrt(J)` scaling.
#' @details
#' Outside the margin a negative Delta_2 chooses information and a positive
#' Delta_2 chooses bootstrap. Within the margin the MLE uses the sign of Delta_3
#' in the same way; zero Delta_3 makes no choice. The WLE makes no choice within
#' the margin because its third-order coefficient is unavailable.
#'
#' The following scope statements concern logistic designs and the MLE
#' conditional on non-extreme responses. The WLE part of the note has
#' not been refereed.
#' \enumerate{
#'   \item On a compact set of abilities containing no zero of Delta_2,
#'   choosing by its sign at the estimate gives a scaled risk within an
#'   exponentially small amount of the smaller of the two fixed scaled risks.
#'   \item Near a switch, any rule switching between the two SEs at the estimate
#'   has the smaller fixed scaled risk plus or minus a term of order `J^(-3/2)`.
#'   Its sign is the sign of c G' Delta_2' at the switch and cannot be chosen by
#'   the user. Always reporting one SE costs at most |D_J|, of order `J^(-2)`.
#'   A margin moves the switch to the margin band's edge and does not remove
#'   this term.
#'   \item For a 2P design with half the items (a, b) = (1, -1/2) and half
#'   (2, 1), at ability 0.01805334808 near a zero of Delta_2, sign selection
#'   exceeds both fixed scaled risks by 18.1%, 16.9%, 9.1% and 4.6% of
#'   the smaller risk at J = 50, 100, 400 and 1600, respectively.
#'   \item There is no guarantee at a fixed test length. Every asymptotic
#'   statement assumes design averages that do not change with J, such as
#'   identical items or a fixed item block repeated.
#' }
#'
#' Exact enumeration stops above 100000 atoms or lattice width 1e6; the package
#' stops instead of approximating. The MLE requires at least two items. A missing
#' or non-unique root causes an error naming T.
#' WLE roots are counted on a grid and at function extrema located by derivative
#' sign changes and searches at sampled local extrema of the derivative.
#' Features narrower than the grid spacing can still be missed; this count
#' does not certify global uniqueness.
#' @template numerical-accuracy
#' @export
#' @examples
#' se_select(rep(1, 4), rep(0, 4), 0, "MLE", margin = 0)
se_select <- function(a, b, theta, estimator, margin) {
  validate_design(a, b, theta)
  validate_law_estimator(a, estimator)
  validate_margin(margin)
  eta <- checked_logits(a * (theta - b))
  info <- information_se(a, eta)
  boot <- bootstrap_se(a, reference_coordinates(b)$b, eta, estimator)
  coefficients <- rule_coefficients(a, eta, estimator)
  d2 <- coefficients$delta2
  d3 <- coefficients$delta3
  decision <- selection_decision(d2, d3, estimator, margin)
  selected <- switch(decision$choice, info = info, boot = boot, none = NA_real_)
  list(choice = decision$choice, se = selected, se_info = info, se_boot = boot,
       delta2 = d2, delta3 = d3, within_margin = decision$within_margin,
       reason = decision$reason)
}
