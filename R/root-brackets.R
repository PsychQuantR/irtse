root_error <- function(T, estimator, detail) {
  stop(sprintf("%s root for T = %.17g: %s.", estimator, T, detail), call. = FALSE)
}

root_brackets <- function(T, a, b, estimator, interval) {
  if (estimator == "MLE") {
    low <- rep(interval[[1L]], length(T))
    high <- rep(interval[[2L]], length(T))
    fl <- estimate_equation(low, a, b, estimator, T)$value
    fh <- estimate_equation(high, a, b, estimator, T)$value
    bad <- which(fl < 0 | fh > 0)
    if (length(bad)) root_error(T[bad[[1L]]], estimator, "no sign change")
    return(list(low = low, high = high))
  }
  grid <- wle_grid(a, b, interval)
  derivatives <- grid_equation(grid, a, b, derivative = TRUE)
  extrema <- which(sign(derivatives[-length(grid)]) * sign(derivatives[-1L]) < 0)
  if (length(extrema)) {
    roots <- refine_roots(rep(0, length(extrema)), grid[extrema],
                          grid[extrema + 1L], a, b, estimator, derivative = TRUE)
    grid <- sort(unique(c(grid, roots)))
  }
  brackets <- vapply(T, function(total) {
    f <- grid_equation(grid, a, b, T = total)
    hits <- which(f == 0)
    changes <- which(sign(f[-length(f)]) * sign(f[-1L]) < 0)
    if (length(hits) + length(changes) != 1L) {
      root_error(total, estimator, "zero or multiple roots on the grid and its extrema")
    }
    if (length(hits)) return(rep(grid[hits], 2L))
    grid[c(changes, changes + 1L)]
  }, numeric(2))
  list(low = brackets[1L, ], high = brackets[2L, ])
}

wle_grid <- function(a, b, interval) {
  size <- max(2001L, 40 * length(a) + 1)
  grid <- seq(interval[[1L]], interval[[2L]], length.out = size)
  # Local grids resolve sharp transitions of highly discriminating items.
  span <- 4 * (log(max(a)) - log(min(a))) + log(4 * length(a)) + 1
  need_local <- which(a * (diff(interval) / (size - 1)) > span / 1000)
  designs <- unique(data.frame(a = a[need_local], b = b[need_local]))
  local <- lapply(seq_len(nrow(designs)), function(i) {
    designs$b[[i]] + seq(-span, span, length.out = 2001L) / designs$a[[i]]
  })
  points <- c(grid, unlist(local, use.names = FALSE))
  sort(unique(points[points >= interval[[1L]] & points <= interval[[2L]]]))
}

grid_equation <- function(grid, a, b, T = 0, derivative = FALSE) {
  # Evaluate every point, with bounded temporary matrix sizes.
  size <- max(1L, floor(1e6 / length(a)))
  starts <- seq(1, length(grid), by = size)
  unlist(lapply(starts, function(start) {
    indices <- seq.int(start, min(length(grid), start + size - 1L))
    equation <- estimate_equation(grid[indices], a, b, "WLE", T)
    if (derivative) equation$derivative else equation$value
  }), use.names = FALSE)
}
