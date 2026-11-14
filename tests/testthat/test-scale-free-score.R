# Requirement: Estimators on the sufficient statistic
# Scenario: Score terms below the double range
test_that("MLE roots use log score terms beyond the probability range", {
  s <- 1000
  expected <- log(2) / (2 * s)
  for (B in c(0.7455, 0.746, 0.76, 5, 500)) {
    for (direction in c(-1, 1)) {
      b <- direction * c(B, -B, -B)
      total <- if (direction == -1) s else 2 * s
      root <- tryCatch(irtse:::estimate_atoms(total, rep(s, 3), b, "MLE"),
                       error = conditionMessage)
      expect_true(is.numeric(root) &&
                    abs(root / (direction * expected) - 1) < 1e-6,
                  info = paste("B =", B, "direction =", direction,
                               "returned", root))
    }
  }
})

# Requirement: Estimators on the sufficient statistic; Numerical range
test_that("Bootstrap SE uses the correct root with underflowed probabilities", {
  s <- 1000
  B <- 0.746
  expected <- (B + log(2) / (2 * s)) * sqrt(2) * exp(-s * B / 2)
  expect_lt(abs(se_boot(rep(s, 3), c(-B, B, B), 0, "MLE") / expected - 1),
            1e-6)
})

# Requirement: Estimators on the sufficient statistic; Numerical range
test_that("A bootstrap SE below the double range stops after locating roots", {
  # (5 + log(2) / 2000) * sqrt(2) * exp(-2500) is positive but underflows.
  expect_error(se_boot(rep(1000, 3), c(-5, 5, 5), 0, "MLE"),
               "representable", fixed = TRUE)
})
