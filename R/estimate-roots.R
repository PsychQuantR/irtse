refine_roots <- function(T, low, high, a, b, estimator, derivative = FALSE) {
  evaluate <- function(x, totals) {
    eq <- estimate_equation(x, a, b, estimator, totals)
    if (derivative) list(value = eq$derivative, derivative = rep(0, length(x))) else eq
  }
  fl <- evaluate(low, T)$value
  fh <- evaluate(high, T)$value
  roots <- low + (high - low) / 2
  roots[fl == 0] <- low[fl == 0]
  roots[fh == 0] <- high[fh == 0]
  active <- which(fl != 0 & fh != 0 & low != high)
  for (iteration in seq_len(2200L)) {
    if (!length(active)) return(roots)
    x <- roots[active]
    eq <- evaluate(x, T[active])
    f <- eq$value
    exact <- f == 0
    left <- sign(f) == sign(fl[active]) & !exact
    right <- !left & !exact
    low[active[left]] <- x[left]
    fl[active[left]] <- f[left]
    high[active[right]] <- x[right]
    fh[active[right]] <- f[right]
    middle <- low[active] + (high[active] - low[active]) / 2
    adjacent <- middle == low[active] | middle == high[active]
    done <- exact | adjacent
    best <- ifelse(abs(fl[active]) <= abs(fh[active]),
                   low[active], high[active])
    roots[active[done]] <- ifelse(exact[done], x[done], best[done])
    next_x <- x - f / eq$derivative
    width <- high[active] - low[active]
    safe <- is.finite(next_x) & next_x > low[active] + width / 10 &
      next_x < high[active] - width / 10
    next_x[!safe] <- middle[!safe]
    roots[active[!done]] <- next_x[!done]
    active <- active[!done]
  }
  root_error(T[active[[1L]]], estimator, "refinement did not converge")
}

estimate_atoms <- function(T, a, b, estimator) {
  L <- (log(4 * length(a)) + log(max(a)) - log(min(a))) / min(a) + 1
  interval <- finite_value(c(min(b) - L, max(b) + L), "root interval")
  totals <- unique(T)
  bracket <- root_brackets(totals, a, b, estimator, interval)
  roots <- refine_roots(totals, bracket$low, bracket$high, a, b, estimator)
  roots[match(T, totals)]
}
