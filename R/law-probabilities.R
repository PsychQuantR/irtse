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

atom_log_probabilities <- function(plan, eta) {
  probabilities <- log_probabilities(eta)
  mass <- if (plan$type == "lattice") {
    lattice_probabilities(plan$weights, probabilities$p, probabilities$q)[plan$keep]
  } else {
    pattern_probabilities(plan$patterns, probabilities$p, probabilities$q)
  }
  total <- log_sum(mass)
  if (!is.finite(total)) representation_error("law of T at this theta")
  mass - total
}

atom_probabilities <- function(plan, b, theta) {
  exp(atom_log_probabilities(plan, plan$a * (theta - b)))
}

law_sd <- function(plan, estimates, eta) {
  weighted_sd(atom_log_probabilities(plan, eta), estimates)
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
