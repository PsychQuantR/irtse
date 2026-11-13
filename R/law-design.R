enumeration_error <- function(reason) {
  stop(paste0("exact enumeration is not available for this design: ", reason,
              "; the package does not approximate."), call. = FALSE)
}

check_atom_count <- function(n, limit, caller) {
  if (!is.finite(n) || n > limit) {
    enumeration_error(paste0(
      "the law of T has ", format(n, scientific = FALSE, big.mark = ""),
      " atoms, above the ", caller, " limit of ",
      format(limit, scientific = FALSE, big.mark = "")
    ))
  }
  invisible(NULL)
}

lattice_weights <- function(a) {
  ratios <- a / min(a)
  for (m in seq_len(1000L)) {
    scaled <- m * ratios
    if (all(is.finite(scaled)) &&
        all(abs(scaled - round(scaled)) <= 8 * .Machine$double.eps * scaled)) {
      return(list(weights = round(scaled), step = min(a) / m))
    }
  }
  NULL
}

reachable_points <- function(weights) {
  reachable <- TRUE
  for (w in weights) {
    old <- seq_along(reachable)
    next_reachable <- rep(FALSE, length(reachable) + w)
    next_reachable[old] <- reachable
    next_reachable[old + w] <- next_reachable[old + w] | reachable
    reachable <- next_reachable
  }
  which(reachable) - 1L
}

atom_plan <- function(a, estimator, limit, caller) {
  truncated <- estimator == "MLE"
  lattice <- lattice_weights(a)
  if (is.null(lattice)) return(pattern_plan(a, truncated, limit, caller))
  width <- sum(lattice$weights)
  if (!is.finite(width) || width > 1e6) {
    enumeration_error(paste0(
      "the lattice width is ",
      format(width, scientific = FALSE, big.mark = ""),
      ", above the limit 1e6 (largest lattice weight ",
      format(max(lattice$weights), scientific = FALSE, big.mark = ""), ")"
    ))
  }
  points <- reachable_points(lattice$weights)
  if (truncated) points <- points[points > 0 & points < width]
  check_atom_count(length(points), limit, caller)
  list(type = "lattice", T = points * lattice$step, keep = points + 1L,
       weights = lattice$weights, truncated = truncated)
}

pattern_plan <- function(a, truncated, limit, caller) {
  count <- 2^length(a)
  check_atom_count(count - 2 * truncated, limit, caller)
  ids <- if (truncated) seq_len(count - 2L) else 0:(count - 1L)
  patterns <- vapply(seq_along(a), function(j) {
    (ids %/% 2^(j - 1L)) %% 2 == 1
  }, logical(length(ids)))
  list(type = "patterns", T = as.vector(patterns %*% a),
       patterns = patterns, truncated = truncated)
}
