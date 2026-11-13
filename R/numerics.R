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

signed_log_mean <- function(log_terms, factors) {
  keep <- factors != 0
  if (!any(keep)) return(c(sign = 0, log_abs = -Inf))
  logs <- log_terms[keep] + log(abs(factors[keep]))
  signs <- sign(factors[keep])
  # Cancel equal opposite terms before exponentiation, so a tiny remaining
  # term is not rounded away between large symmetric contributions.
  if (any(signs > 0) && any(signs < 0)) {
    unique_logs <- unique(logs)
    signs <- as.vector(rowsum(signs, match(logs, unique_logs), reorder = FALSE))
    logs <- unique_logs + log(abs(signs))
    keep <- signs != 0
    logs <- logs[keep]
    signs <- sign(signs[keep])
  }
  if (!length(logs)) return(c(sign = 0, log_abs = -Inf))
  largest <- max(logs)
  if (!is.finite(largest)) representation_error("design average")
  total <- sum(signs * exp(logs - largest))
  c(sign = sign(total), log_abs = largest + log(abs(total)) - log(length(factors)))
}

# One source for information, cumulants and discrimination-scaled k_r / k_2.
# Scaling a by exp(log_scale) keeps the coefficient polynomials bounded.
# The final powers of a and k2 are restored together on the log scale.
design_ratios <- function(a, b, theta) {
  eta <- a * (theta - b)
  log_u <- plogis(eta, log.p = TRUE) + plogis(-eta, log.p = TRUE)
  u <- exp(log_u)
  s <- plogis(-eta) - plogis(eta)
  factors <- list(rep(1, length(a)), s, 1 - 6 * u, s * (1 - 12 * u),
                  1 - 30 * u + 120 * u^2, s * (1 - 60 * u + 360 * u^2))
  log_a <- log(a)
  log_weights <- 2 * log_a + log_u
  largest <- max(log_weights)
  if (!is.finite(largest)) representation_error("information")
  relative <- log_weights - largest
  log_total <- log_sum(relative)
  # The fifth weighted moment of a bounds every scaled ratio through k7.
  # A completely saturated item therefore cannot set the working scale.
  log_scale <- (log_sum(relative + 5 * log_a) - log_total) / 5
  averages <- vapply(2:7, function(r) {
    signed_log_mean(r * log_a + log_u, factors[[r - 1L]])
  }, numeric(2))
  log_k2 <- averages["log_abs", 1L]
  if (!is.finite(log_k2)) representation_error("information")
  ratios <- vapply(2:7, function(r) {
    average <- signed_log_mean(relative + (r - 2) * (log_a - log_scale),
                               factors[[r - 1L]])
    exp_value(average[["log_abs"]] - log_total + log(length(a)),
              average[["sign"]], "cumulant ratio")
  }, numeric(1))
  names(ratios) <- paste0("k", 2:7)
  list(ratios = finite_value(ratios, "cumulant ratios"), log_k2 = log_k2,
       log_scale = log_scale, averages = averages)
}
