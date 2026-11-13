enumeration_error <- function(limit) {
  stop(paste0("exact enumeration is not available: limit ", limit,
              " exceeded; the package does not approximate."), call. = FALSE)
}

check_atom_count <- function(n, limit) {
  if (!is.finite(n) || n > limit) {
    enumeration_error(format(limit, scientific = FALSE, trim = TRUE))
  }
  invisible(NULL)
}

lattice_weights <- function(a) {
  ratios <- a / min(a)
  for (m in seq_len(1000L)) {
    scaled <- m * ratios
    if (all(is.finite(scaled)) &&
        all(abs(scaled - round(scaled)) <= 1e-12 * scaled)) {
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

atom_plan <- function(a, estimator, limit) {
  truncated <- estimator == "MLE"
  if (all(a == a[[1L]])) {
    width <- length(a)
    if (width > 1e6) enumeration_error("1e6 (1000000)")
    check_atom_count(width + 1 - 2 * truncated, limit)
    points <- if (truncated) seq_len(width - 1L) else 0:width
    return(list(type = "1P", T = a[[1L]] * points,
                keep = points + 1L, truncated = truncated))
  }
  lattice <- lattice_weights(a)
  if (is.null(lattice)) return(pattern_plan(a, truncated, limit))
  width <- sum(lattice$weights)
  if (!is.finite(width) || width > 1e6) {
    enumeration_error("1e6 (1000000)")
  }
  points <- reachable_points(lattice$weights)
  if (truncated) points <- points[points > 0 & points < width]
  check_atom_count(length(points), limit)
  list(type = "lattice", T = points * lattice$step, keep = points + 1L,
       weights = lattice$weights, truncated = truncated)
}

pattern_plan <- function(a, truncated, limit) {
  count <- 2^length(a)
  check_atom_count(count - 2 * truncated, limit)
  ids <- if (truncated) seq_len(count - 2L) else 0:(count - 1L)
  patterns <- vapply(seq_along(a), function(j) {
    (ids %/% 2^(j - 1L)) %% 2 == 1
  }, logical(length(ids)))
  list(type = "patterns", T = as.vector(patterns %*% a),
       patterns = patterns, truncated = truncated)
}
