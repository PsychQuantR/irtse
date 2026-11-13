estimate_equation <- function(theta, a, b, estimator) {
  eta <- outer(theta, a) - rep(a * b, each = length(theta))
  p <- plogis(eta)
  u <- p * plogis(-eta)
  information <- as.vector(u %*% a^2)
  score <- -as.vector(p %*% a)
  derivative <- -information
  if (estimator == "WLE") {
    third <- as.vector((u * (1 - 2 * p)) %*% a^3)
    fourth <- as.vector((u * (1 - 6 * u)) %*% a^4)
    score <- score + third / (2 * information)
    derivative <- derivative + (fourth / information -
                                  (third / information)^2) / 2
  }
  if (any(!is.finite(score)) || any(!is.finite(derivative))) {
    stop("a and b yield a non-finite estimating equation.", call. = FALSE)
  }
  list(value = score, derivative = derivative)
}

root_error <- function(T, estimator, detail) {
  stop(sprintf("%s root for T = %.17g: %s.", estimator, T, detail),
       call. = FALSE)
}

root_brackets <- function(T, a, b, estimator, tau) {
  if (estimator == "MLE") {
    low <- rep(-tau, length(T))
    high <- rep(tau, length(T))
    fl <- T + estimate_equation(-tau, a, b, estimator)$value
    fh <- T + estimate_equation(tau, a, b, estimator)$value
    bad <- which(fl < 0 | fh > 0)
    if (length(bad)) root_error(T[bad[[1L]]], estimator, "no sign change")
    return(list(low = low, high = high))
  }
  grid <- wle_grid(a, b, tau)
  g <- grid_equation(grid, a, b)
  brackets <- vapply(T, function(total) {
    f <- total + g
    hits <- which(f == 0)
    changes <- which(sign(f[-length(f)]) * sign(f[-1L]) < 0)
    if (length(hits) + length(changes) != 1L) {
      root_error(total, estimator, "zero or multiple roots on the grid")
    }
    if (length(hits)) return(rep(grid[hits], 2L))
    grid[c(changes, changes + 1L)]
  }, numeric(2))
  list(low = brackets[1L, ], high = brackets[2L, ])
}

wle_grid <- function(a, b, tau) {
  size <- max(2001L, 40 * length(a) + 1)
  grid <- seq(-tau, tau, length.out = size)
  # Local grids resolve sharp transitions of highly discriminating items.
  # The span also covers where their higher derivatives dominate other items.
  span <- 4 * log(max(a) / min(a)) + log(4 * length(a)) + 1
  need_local <- which(a * (2 * tau / (size - 1)) > span / 1000)
  designs <- unique(data.frame(a = a[need_local], b = b[need_local]))
  local <- lapply(seq_len(nrow(designs)), function(i) {
    designs$b[[i]] + seq(-span, span, length.out = 2001L) / designs$a[[i]]
  })
  points <- c(grid, unlist(local, use.names = FALSE))
  sort(unique(points[points >= -tau & points <= tau]))
}

grid_equation <- function(grid, a, b) {
  # Evaluate every point, with bounded temporary matrix sizes.
  size <- max(1L, floor(1e6 / length(a)))
  starts <- seq(1, length(grid), by = size)
  unlist(lapply(starts, function(start) {
    indices <- seq.int(start, min(length(grid), start + size - 1L))
    estimate_equation(grid[indices], a, b, "WLE")$value
  }), use.names = FALSE)
}
