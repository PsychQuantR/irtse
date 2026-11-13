estimate_equation <- function(theta, a, b, estimator, T = 0) {
  eta <- sweep(outer(theta, b, "-"), 2L, a, "*")
  p <- plogis(eta)
  q <- plogis(-eta)
  positive <- eta > 0
  weighted <- function(x, powers) rowSums(sweep(x, 2L, powers, "*"))
  score <- (T - weighted(positive, a)) + weighted(ifelse(positive, q, 0), a) -
    weighted(ifelse(positive, 0, p), a)
  log_weights <- sweep(plogis(eta, log.p = TRUE) +
                         plogis(-eta, log.p = TRUE), 2L, 2 * log(a), "+")
  largest <- apply(log_weights, 1L, max)
  if (any(!is.finite(largest))) representation_error("estimating equation")
  weights <- exp(log_weights - largest)
  total <- rowSums(weights)
  information <- exp(largest + log(total))
  if (estimator == "MLE" && any(score == 0 & information == 0)) {
    representation_error("estimating equation tails")
  }
  derivative <- -information
  if (estimator == "WLE") {
    third <- weighted(weights * (q - p), a) / total
    fourth <- weighted(weights * (1 - 6 * p * q), a^2) / total
    score <- score + third / 2
    derivative <- derivative + (fourth - third^2) / 2
  }
  list(value = finite_value(score, "estimating equation"),
       derivative = finite_value(derivative, "estimating derivative"))
}
