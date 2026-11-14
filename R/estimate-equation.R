estimate_equation <- function(theta, a, b, estimator, T = 0) {
  eta <- sweep(outer(theta, b, "-"), 2L, a, "*")
  probabilities <- log_probabilities(eta)
  terms <- score_terms(eta, a, T, probabilities)
  log_weights <- sweep(probabilities$p + probabilities$q, 2L, 2 * log(a), "+")
  largest <- apply(log_weights, 1L, max)
  if (any(!is.finite(largest))) representation_error("estimating equation")
  weights <- exp(log_weights - largest)
  total <- rowSums(weights)
  log_information <- largest + log(total)
  correction <- correction_derivative <- numeric(length(theta))
  if (estimator == "WLE") {
    adjustment <- wle_adjustment(weights, total, eta, a, probabilities)
    correction <- adjustment$value
    correction_derivative <- adjustment$derivative
  }
  logs <- cbind(terms$logs, log(abs(correction)))
  signs <- cbind(terms$signs, sign(correction))
  scale <- apply(logs, 1L, max)
  scale[scale == -Inf] <- 0 # An exactly zero function needs no rescaling.
  value <- compensated_rows(signs * exp(logs - scale))
  scaled_derivative <- -exp(log_information - scale) +
    sign(correction_derivative) * exp(log(abs(correction_derivative)) - scale)
  list(value = finite_value(sign(value) * exp(log(abs(value)) + scale),
                            "estimating equation"),
       derivative = finite_value(-exp(log_information) + correction_derivative,
                                 "estimating derivative"),
       scaled_value = value, scaled_derivative = scaled_derivative,
       log_scale = scale)
}

score_terms <- function(eta, a, T, probabilities) {
  positive <- eta > 0
  central <- abs(eta) < 1
  constant <- -sweep(positive, 2L, a, "*")
  halves <- matrix(rep(a / 2, each = nrow(eta)), nrow = nrow(eta))
  constant[central] <- -halves[central]
  deviation <- matrix(0, nrow(eta), ncol(eta))
  # Keep the deviation from 1/2 separately even when probabilities round to it.
  deviation[central] <- -halves[central] * tanh(eta[central] / 2)
  residual <- compensated_rows(cbind(constant, deviation), T)
  logs <- sweep(ifelse(positive, probabilities$q, probabilities$p), 2L, log(a), "+")
  logs[central] <- -Inf
  list(logs = cbind(log(abs(residual)), logs),
       signs = cbind(sign(residual), ifelse(positive, 1, -1)))
}

wle_adjustment <- function(weights, total, eta, a, probabilities) {
  weighted <- function(x, powers) rowSums(sweep(x, 2L, powers, "*"))
  third <- weighted(weights * (-tanh(eta / 2)), a) / total
  u <- exp(probabilities$p) * exp(probabilities$q)
  fourth <- weighted(weights * (1 - 6 * u), a^2) / total
  list(value = third / 2, derivative = (fourth - third^2) / 2)
}

compensated_rows <- function(terms, initial = 0) {
  total <- rep_len(initial, nrow(terms))
  compensation <- numeric(nrow(terms))
  for (j in seq_len(ncol(terms))) {
    term <- terms[, j]
    next_total <- total + term
    larger <- abs(total) >= abs(term)
    error <- ifelse(larger, (total - next_total) + term,
                    (term - next_total) + total)
    compensation <- compensation + error
    total <- next_total
  }
  total + compensation
}
