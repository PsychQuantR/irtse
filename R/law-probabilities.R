lattice_probabilities <- function(weights, log_p, log_q) {
  mass <- 0
  for (j in seq_along(weights)) {
    indices <- seq_along(mass)
    next_mass <- rep(-Inf, length(mass) + weights[[j]])
    next_mass[indices] <- mass + log_q[[j]]
    shifted <- indices + weights[[j]]
    next_mass[shifted] <- log_add(next_mass[shifted], mass + log_p[[j]])
    mass <- next_mass
  }
  mass
}

pattern_probabilities <- function(patterns, log_p, log_q) {
  mass <- numeric(nrow(patterns))
  for (j in seq_along(log_p)) {
    mass <- mass + ifelse(patterns[, j], log_p[[j]], log_q[[j]])
  }
  mass
}

atom_log_probabilities <- function(plan, a, b, theta) {
  eta <- a * (theta - b)
  log_p <- plogis(eta, log.p = TRUE)
  log_q <- plogis(-eta, log.p = TRUE)
  mass <- if (plan$type == "lattice") {
    lattice_probabilities(plan$weights, log_p, log_q)[plan$keep]
  } else {
    pattern_probabilities(plan$patterns, log_p, log_q)
  }
  total <- log_sum(mass)
  if (!is.finite(total)) representation_error("law of T at this theta")
  mass - total
}

atom_probabilities <- function(plan, a, b, theta) {
  exp(atom_log_probabilities(plan, a, b, theta))
}

law_sd <- function(plan, estimates, a, b, theta) {
  weighted_sd(atom_log_probabilities(plan, a, b, theta), estimates)
}

weighted_sd <- function(logs, estimates) {
  # Centre relative to the modal estimate, preserving translation invariance
  # when all estimates are far from zero.
  offsets <- estimates - estimates[[which.max(logs)]]
  centre <- sum(exp(logs) * offsets)
  distances <- abs(offsets - centre)
  if (all(distances == 0)) return(0)
  variance <- log_sum(logs + 2 * log(distances))
  exp_value(variance / 2, what = "bootstrap SE")
}
