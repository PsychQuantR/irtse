# Requirement: Estimators on the sufficient statistic
test_that("MLE roots for identical Rasch items", {
  h0 <- sqrt((36 * log(5)^2 + 90 * log(2)^2) / 31)

  expect_equal(
    se_boot(a = rep(1, 6), b = rep(0, 6), theta = 0, estimator = "MLE"),
    h0 / sqrt(6),
    tolerance = 1e-12
  )
})

# Requirement: Estimators on the sufficient statistic
test_that("WLE roots for identical Rasch items", {
  totals <- 0:4
  roots <- qlogis((totals + 1/2) / 5)
  weights <- choose(4, totals) / 16
  root_mean <- sum(weights * roots)
  root_sd <- sqrt(sum(weights * (roots - root_mean)^2))

  expect_equal(
    se_boot(a = rep(1, 4), b = rep(0, 4), theta = 0, estimator = "WLE"),
    root_sd,
    tolerance = 1e-12
  )
})

# Requirement: Exact law of the sufficient statistic
test_that("1P with a common discrimination other than one", {
  expect_equal(
    se_boot(
      a = rep(2, 5), b = c(-1, 0, 0.5, 1, 2), theta = 0.3,
      estimator = "MLE"
    ),
    se_boot(
      a = rep(1, 5), b = 2 * c(-1, 0, 0.5, 1, 2), theta = 0.6,
      estimator = "MLE"
    ) / 2,
    tolerance = 1e-10
  )
})

# Requirement: Exact law of the sufficient statistic
test_that("Lattice path and pattern path agree", {
  a <- c(1, 2, 3, 1.5)
  pattern_a <- c(1, 2, 3, 1.5 * (1 + 1e-9))
  b <- c(-0.5, 0, 0.4, 1)

  for (estimator in c("MLE", "WLE")) {
    expect_equal(
      se_boot(a = a, b = b, theta = 0.2, estimator = estimator),
      se_boot(a = pattern_a, b = b, theta = 0.2, estimator = estimator),
      tolerance = 1e-6
    )
  }
})

# Requirement: Information-based standard error
test_that("Six identical Rasch items at their difficulty", {
  expect_equal(
    se_info(a = rep(1, 6), b = rep(0, 6), theta = 0),
    sqrt(2/3),
    tolerance = 1e-12
  )
})

# Requirement: Exact parametric-bootstrap standard error
test_that("Bootstrap SE is the scaled target divided by sqrt(J)", {
  expect_equal(
    se_boot(a = rep(1, 6), b = rep(0, 6), theta = 0, estimator = "MLE"),
    risk_diff(
      a = rep(1, 6), b = rep(0, 6), theta = 0, estimator = "MLE"
    )$h_true / sqrt(6),
    tolerance = 1e-12
  )
})

# Requirement: Limits of exact enumeration
test_that("Pattern path too large for se_boot", {
  a <- 1 + sqrt(2) * (1:18) / 100

  # Both estimators must reject the design before enumerating its patterns.
  for (estimator in c("MLE", "WLE")) {
    expect_error(
      se_boot(a = a, b = rep(0, 18), theta = 0, estimator = estimator),
      "100000"
    )
    expect_error(
      se_boot(a = a, b = rep(0, 18), theta = 0, estimator = estimator),
      "(?i)exact enumeration", perl = TRUE
    )
  }
})

# Requirement: Limits of exact enumeration
test_that("Pattern path too large for risk_diff", {
  a <- 1 + sqrt(2) * (1:12) / 100

  for (estimator in c("MLE", "WLE")) {
    expect_error(
      risk_diff(a = a, b = rep(0, 12), theta = 0, estimator = estimator),
      "2000"
    )
    expect_error(
      risk_diff(a = a, b = rep(0, 12), theta = 0, estimator = estimator),
      "(?i)exact enumeration", perl = TRUE
    )
  }
})

# Requirement: Input validation at the boundary
test_that("Invalid inputs", {
  expect_error(se_info(a = c(1, -1), b = c(0, 0), theta = 0), "a")
  expect_error(se_info(a = c(1, 1), b = c(0, NA), theta = 0), "b")
  expect_error(se_info(a = c(1, 1), b = 0, theta = 0), "b")
  expect_error(
    se_info(a = c(1, 1), b = c(0, 0), theta = c(0, 1)), "theta"
  )
  expect_error(
    se_boot(a = c(1, 1), b = c(0, 0), theta = 0), "estimator"
  )
  expect_error(
    se_boot(a = c(1, 1), b = c(0, 0), theta = 0, estimator = "EAP"),
    "estimator"
  )
  expect_error(se_boot(a = 1, b = 0, theta = 0, estimator = "MLE"), "MLE")
})

# Requirement: Package metadata and documentation
test_that("README statements", {
  readme_path <- testthat::test_path("..", "..", "README.md")
  if (!file.exists(readme_path)) {
    skip("README not available in installed package")
  }

  # The README describes a public API; require that API before checking prose.
  for (name in c(
    "se_info", "se_boot", "design_cumulants", "delta2", "delta3",
    "risk_diff", "se_select"
  )) {
    getExportedValue("irtse", name)
  }
  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")

  expect_match(readme, "unreleased", fixed = TRUE)
  expect_match(readme, "working name", fixed = TRUE)
  expect_match(readme, "catR::semTheta", fixed = TRUE)
})
