# The number of abilities whose estimating equation is evaluated together: at
# most 1e6 / J, and at least one. A block of n abilities makes temporary matrices
# of n rows and at most 2J + 2 columns (the residual's terms in score_terms()
# have 2J, the scaled terms of the equation J + 2), so their size does not
# depend on the number of abilities.
equation_block <- function(J) max(1L, as.integer(floor(1e6 / J)))

# The estimating equation at each ability `theta` with its total `T` (recycled
# to the length of `theta`). Every row is computed from that row alone, so the
# rows are evaluated in blocks of equation_block() rows, which bounds the memory
# by the block, not by the number of abilities. The three checks run on all rows
# after the last block, in the order of an evaluation of all rows at once, so
# the blocks change neither a returned double nor the error a call stops with.
estimate_equation <- function(theta, a, b, estimator, T = 0) {
  n <- length(theta)
  T <- rep_len(T, n)
  size <- equation_block(length(a))
  starts <- if (n == 0L) 1L else seq.int(1L, n, by = size)
  parts <- lapply(starts, function(start) {
    rows <- seq_len(n)[seq.int(start, length.out = min(size, n - start + 1L))]
    equation_rows(theta[rows], a, b, estimator, T[rows])
  })
  part <- function(name) unlist(lapply(parts, `[[`, name), use.names = FALSE)
  if (any(!is.finite(part("largest")))) representation_error("estimating equation")
  value <- part("value")
  scale <- part("scale")
  list(value = finite_value(sign(value) * exp(log(abs(value)) + scale),
                            "estimating equation"),
       derivative = finite_value(-exp(part("log_information")) +
                                   part("correction_derivative"),
                                 "estimating derivative"),
       scaled_value = value, scaled_derivative = part("scaled_derivative"),
       log_scale = scale)
}

# The quantities of the estimating equation for one block of rows, before the
# checks of estimate_equation(). A row whose largest log weight is not finite
# gives values that are not finite here; the caller stops on it.
equation_rows <- function(theta, a, b, estimator, T) {
  eta <- sweep(outer(theta, b, "-"), 2L, a, "*")
  probabilities <- log_probabilities(eta)
  terms <- score_terms(eta, a, T, probabilities)
  log_weights <- sweep(probabilities$p + probabilities$q, 2L, 2 * log(a), "+")
  largest <- apply(log_weights, 1L, max)
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
  list(largest = largest, log_information = log_information, value = value,
       scale = scale, scaled_derivative = scaled_derivative,
       correction_derivative = correction_derivative)
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
