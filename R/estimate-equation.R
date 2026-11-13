estimate_equation <- function(theta, a, b, estimator, T = 0) {
  eta <- sweep(outer(theta, b, "-"), 2L, a, "*")
  probabilities <- log_probabilities(eta)
  p <- plogis(eta)
  q <- plogis(-eta)
  positive <- eta > 0
  weighted <- function(x, powers) rowSums(sweep(x, 2L, powers, "*"))
  score <- compensated_score(eta, a, T, p, q)
  log_weights <- sweep(probabilities$p + probabilities$q, 2L, 2 * log(a), "+")
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
    third <- weighted(weights * (-tanh(eta / 2)), a) / total
    fourth <- weighted(weights * (1 - 6 * p * q), a^2) / total
    score <- score + third / 2
    derivative <- derivative + (fourth - third^2) / 2
  }
  list(value = finite_value(score, "estimating equation"),
       derivative = finite_value(derivative, "estimating derivative"))
}

compensated_score <- function(eta, a, T, p, q) {
  positive <- eta > 0
  central <- abs(eta) < 1
  constant <- -sweep(positive, 2L, a, "*")
  tail <- sweep(ifelse(positive, q, -p), 2L, a, "*")
  # Keep the deviation from 1/2 separately: plogis alone rounds it away
  # for small logits, even when compensated accumulation retains the tails.
  halves <- matrix(rep(a / 2, each = nrow(eta)), nrow = nrow(eta))
  constant[central] <- -halves[central]
  tail[central] <- -halves[central] * tanh(eta[central] / 2)
  terms <- cbind(constant, tail)
  total <- rep_len(T, nrow(eta))
  compensation <- numeric(nrow(eta))
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
