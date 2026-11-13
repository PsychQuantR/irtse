lattice_probabilities <- function(weights, p, q) {
  mass <- 1
  for (j in seq_along(weights)) {
    indices <- seq_along(mass)
    next_mass <- numeric(length(mass) + weights[[j]])
    next_mass[indices] <- mass * q[[j]]
    shifted <- indices + weights[[j]]
    next_mass[shifted] <- next_mass[shifted] + mass * p[[j]]
    mass <- next_mass
  }
  mass
}

pattern_probabilities <- function(patterns, p, q) {
  mass <- rep(1, nrow(patterns))
  for (j in seq_along(p)) {
    mass <- mass * ifelse(patterns[, j], p[[j]], q[[j]])
  }
  mass
}

atom_probabilities <- function(plan, a, b, theta) {
  if (plan$type == "1P") {
    mass <- irtcore::lord_wingersky(a[[1L]] * theta, a[[1L]] * b)
  } else {
    eta <- a * (theta - b)
    p <- plogis(eta)
    q <- plogis(-eta)
    mass <- if (plan$type == "lattice") {
      lattice_probabilities(plan$weights, p, q)
    } else {
      pattern_probabilities(plan$patterns, p, q)
    }
  }
  if (plan$type != "patterns") mass <- mass[plan$keep]
  total <- sum(mass)
  if (!is.finite(total) || total <= 0) {
    stop("theta yields a zero or non-finite normalising mass.", call. = FALSE)
  }
  if (plan$truncated) mass <- mass / total
  mass
}

law_sd <- function(plan, estimates, a, b, theta) {
  weights <- atom_probabilities(plan, a, b, theta)
  centre <- sum(weights * estimates)
  sqrt(sum(weights * (estimates - centre)^2))
}
