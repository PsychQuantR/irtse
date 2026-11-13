#' Design averages of logistic cumulants
#'
#' @param a Numeric vector of finite positive item discriminations.
#' @param b Finite numeric difficulties, of the same length as `a`.
#' @param theta Single finite ability, or ability estimate.
#' @return A numeric vector named `k2`, ..., `k7`, with element r equal to
#'   the mean of `a^r` times the rth Bernoulli cumulant. A single item is valid.
#' @template numerical-accuracy
#' @export
#' @examples
#' design_cumulants(1, 0, 0)
design_cumulants <- function(a, b, theta) {
  validate_design(a, b, theta)
  design <- design_ratios(a, b, theta)
  values <- vapply(seq_len(6L), function(i) {
    exp_value(design$averages["log_abs", i], design$averages["sign", i],
              "design cumulant")
  }, numeric(1))
  stats::setNames(values, paste0("k", 2:7))
}
