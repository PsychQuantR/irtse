representation_error <- function(what = "value") {
  stop(paste(what, "is not representable in double precision."), call. = FALSE)
}

finite_value <- function(x, what = "value") {
  if (any(!is.finite(x))) representation_error(what)
  x
}

exp_value <- function(log_abs, sign = 1, what = "value") {
  finite_value(sign, what)
  if (sign == 0) return(0)
  value <- sign * exp(log_abs)
  if (!is.finite(value) || value == 0) representation_error(what)
  unname(value)
}

log_add <- function(x, y) {
  largest <- pmax(x, y)
  result <- largest + log1p(exp(-abs(x - y)))
  result[largest == -Inf] <- -Inf
  result
}

log_sum <- function(x) {
  largest <- max(x)
  if (!is.finite(largest)) return(largest)
  largest + log(sum(exp(x - largest)))
}

unresolved_error <- function(what) {
  stop(paste(what, "is not resolved in double precision."), call. = FALSE)
}

signed_log_mean <- function(log_terms, factors, what = "coefficient polynomial",
                            zero_if_unresolved = FALSE) {
  keep <- factors != 0
  if (!any(keep)) return(c(sign = 0, log_abs = -Inf))
  logs <- log_terms[keep] + log(abs(factors[keep]))
  largest <- max(logs)
  if (!is.finite(largest)) representation_error(what)
  total <- sum(sign(factors[keep]) * exp(logs - largest))
  if (total == 0) {
    if (zero_if_unresolved) return(c(sign = 0, log_abs = -Inf))
    unresolved_error(what)
  }
  c(sign = sign(total), log_abs = largest + log(abs(total)) - log(length(factors)))
}

checked_logits <- function(eta) {
  finite_value(eta, "a logit")
  if (any(eta != 0 & tanh(eta / 2) == 0)) representation_error("a logit")
  eta
}

log_probabilities <- function(eta) {
  eta <- checked_logits(eta)
  list(p = plogis(eta, log.p = TRUE), q = plogis(-eta, log.p = TRUE))
}

log_information <- function(a, eta) {
  probabilities <- log_probabilities(eta)
  log_u <- probabilities$p + probabilities$q
  log_weights <- 2 * log(a) + log_u
  log_k2 <- log_sum(log_weights) - log(length(a))
  if (!is.finite(log_k2)) representation_error("information")
  list(eta = eta, log_u = log_u, log_weights = log_weights, log_k2 = log_k2)
}

odd_unpaired <- function(a, eta) {
  keep <- eta != 0
  if (!any(eta > 0) || !any(eta < 0)) return(keep)
  ordered <- order(a, abs(eta))
  left <- ordered[-length(ordered)]
  right <- ordered[-1L]
  groups <- cumsum(c(TRUE, a[left] != a[right] | abs(eta[left]) != abs(eta[right])))
  for (indices in split(ordered, groups)) {
    positive <- indices[eta[indices] > 0]
    negative <- indices[eta[indices] < 0]
    count <- min(length(positive), length(negative))
    if (count) keep[c(positive[seq_len(count)], negative[seq_len(count)])] <- FALSE
  }
  keep
}

scaled_design_mean <- function(r, log_a, log_u, factors, keep,
                               zero_if_unresolved = FALSE) {
  n <- length(factors)
  keep <- keep & factors != 0
  if (!any(keep)) return(c(sign = 0, log_abs = -Inf))
  log_a <- log_a[keep]
  log_u <- log_u[keep]
  factors <- factors[keep]
  item_logs <- log_u + log(abs(factors))
  anchor <- which.max(r * log_a + item_logs)
  relative_a <- r * (log_a - log_a[[anchor]])
  exponent <- relative_a - item_logs[[anchor]]
  # Multiply u and f before scaling; a common large log(a) must not erase
  # differences between these item terms. Fall back to split log differences
  # only where the separate factors cannot be represented.
  uf <- exp(log_u) * factors
  safe <- abs(uf) >= .Machine$double.xmin &
    exponent < log(.Machine$double.xmax)
  terms <- numeric(length(factors))
  terms[safe] <- exp(exponent[safe]) * uf[safe]
  terms[!safe] <- sign(factors[!safe]) * exp(relative_a[!safe] +
    (log_u[!safe] - log_u[[anchor]]) +
    (log(abs(factors[!safe])) - log(abs(factors[[anchor]]))))
  total <- sum(terms)
  if (total == 0) {
    if (zero_if_unresolved) return(c(sign = 0, log_abs = -Inf))
    unresolved_error(paste0("design average k", r))
  }
  c(sign = sign(total), log_abs = r * log_a[[anchor]] +
      item_logs[[anchor]] + log(abs(total)) - log(n))
}

# Information has one implementation, independent of the higher cumulants.
# Ratios are scaled by a weighted discrimination to keep polynomials bounded.
# zero_if_unresolved is for the zero search of se_diagnostic() only: there an
# average or factor that is not resolved at an ability is a zero of it.
design_ratios <- function(a, eta, zero_if_unresolved = FALSE) {
  information <- log_information(a, eta)
  log_u <- information$log_u
  u <- exp(log_u)
  s <- -tanh(information$eta / 2)
  factors <- list(s, 1 - 6 * u, s * (1 - 12 * u),
                  1 - 30 * u + 120 * u^2, s * (1 - 60 * u + 360 * u^2))
  log_a <- log(a)
  relative <- information$log_weights - max(information$log_weights)
  log_scale <- (log_sum(relative + 5 * log_a) - log_sum(relative)) / 5
  unpaired <- odd_unpaired(a, information$eta)
  averages <- cbind(c(sign = 1, log_abs = information$log_k2),
    vapply(3:7, function(r) {
      scaled_design_mean(r, log_a, log_u, factors[[r - 2L]],
                         if (r %% 2) unpaired else rep(TRUE, length(a)),
                         zero_if_unresolved)
    }, numeric(2)))
  ratios <- vapply(2:7, function(r) {
    exp_value(averages["log_abs", r - 1L] - information$log_k2 -
                (r - 2) * log_scale, averages["sign", r - 1L], "cumulant ratio")
  }, numeric(1))
  names(ratios) <- paste0("k", 2:7)
  list(ratios = finite_value(ratios, "cumulant ratios"),
       log_k2 = information$log_k2, log_scale = log_scale, averages = averages)
}

reference_coordinates <- function(b) {
  reference <- finite_value(min(b) + (max(b) - min(b)) / 2, "reference point")
  list(b = finite_value(b - reference, "shifted difficulties"),
       reference = reference)
}
