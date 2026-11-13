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

# Requirement: Estimators on the sufficient statistic
test_that("WLE equation without a unique root", {
  # An independent 200001-point grid found three sign changes for T = 0
  # and T = 31; the other attainable totals, 1 and 30, each have one.
  expect_error(se_boot(c(1, 30), c(0, 0), 0, "WLE"), "multiple roots")
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
test_that("Equal discriminations at extreme abilities", {
  # Six identical Rasch items: the estimates have closed forms
  # (MLE log(T / (6 - T)), T = 1..5; WLE qlogis((T + 1/2) / 7), T = 0..6) and
  # the law of T is binomial. The weights are formed on the log scale with
  # log P = -log1p(exp(-theta)) and log Q = -theta - log1p(exp(-theta)), so
  # this reference stays exact at theta = 30 and 40. Relative differences are
  # checked explicitly because expect_equal() compares absolutely when the
  # values are below its tolerance.
  reference_sd <- function(theta, totals, estimates, J = 6) {
    log_p <- -log1p(exp(-theta))
    log_q <- -theta - log1p(exp(-theta))
    log_w <- lchoose(J, totals) + totals * log_p + (J - totals) * log_q
    w <- exp(log_w - max(log_w))
    w <- w / sum(w)
    centre <- sum(w * estimates)
    sqrt(sum(w * (estimates - centre)^2))
  }
  for (theta in c(30, 40)) {
    mle <- reference_sd(theta, 1:5, log((1:5) / (6 - (1:5))))
    wle <- reference_sd(theta, 0:6, qlogis(((0:6) + 1/2) / 7))
    x_mle <- se_boot(rep(1, 6), rep(0, 6), theta, "MLE")
    x_wle <- se_boot(rep(1, 6), rep(0, 6), theta, "WLE")

    expect_true(is.finite(x_mle) && x_mle > 0)
    expect_true(is.finite(x_wle) && x_wle > 0)
    expect_lt(abs(x_mle / mle - 1), 1e-8)
    expect_lt(abs(x_wle / wle - 1), 1e-8)
  }
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

  # The message must identify the atom count as the reason for rejection.
  expect_error(
    se_boot(a = a, b = rep(0, 18), theta = 0, estimator = "MLE"), "atom"
  )
})

# Requirement: Limits of exact enumeration
test_that("Lattice limits", {
  # An unmatched error must not prevent the remaining limit checks from running.
  check_error_message <- function(code) {
    tryCatch(force(code), simpleError = function(error) {
      fail(conditionMessage(error))
    })
  }

  check_error_message(expect_error(
    se_boot(c(1, 1000001/1000), c(0, 0), 0, "MLE"),
    "lattice width", fixed = TRUE
  ))
  check_error_message(expect_error(
    se_boot(c(1, 1000001/1000), c(0, 0), 0, "MLE"),
    "1001001", fixed = TRUE
  ))
  check_error_message(expect_error(
    se_boot(c(1, 1000001/1000), c(0, 0), 0, "MLE"),
    "1e6", fixed = TRUE
  ))
  check_error_message(expect_error(
    se_boot(rep(c(1, 997), 600), rep(0, 1200), 0, "WLE"),
    "atom", fixed = TRUE
  ))
  check_error_message(expect_error(
    se_boot(rep(c(1, 997), 600), rep(0, 1200), 0, "WLE"),
    "100000", fixed = TRUE
  ))
  expect_equal(
    se_boot(c(1, 999), c(0, 0), 0, "MLE"), 0.007603208,
    tolerance = 1e-6
  )
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

  # The message must identify the atom count as the reason for rejection.
  expect_error(
    risk_diff(a = a, b = rep(0, 12), theta = 0, estimator = "MLE"), "atom"
  )
})

# Requirement: Input validation at the boundary
test_that("Invalid inputs", {
  expect_error(se_info(a = c(1, -1), b = c(0, 0), theta = 0), "^a must")
  expect_error(se_info(a = c(1, 1), b = c(0, NA), theta = 0), "^b must")
  expect_error(se_info(a = c(1, 1), b = 0, theta = 0), "^b must")
  expect_error(
    se_info(a = c(1, 1), b = c(0, 0), theta = c(0, 1)), "^theta must"
  )
  expect_error(
    se_boot(a = c(1, 1), b = c(0, 0), theta = 0), "^estimator must"
  )
  expect_error(
    se_boot(a = c(1, 1), b = c(0, 0), theta = 0, estimator = "EAP"),
    "^estimator must"
  )
  expect_error(se_boot(a = 1, b = 0, theta = 0, estimator = "MLE"), "MLE")
})

# Requirement: Input validation at the boundary
test_that("Validation in every exported function", {
  expect_error(se_info(a = c(1, -1), b = c(0, 0), theta = 0), "^a must")
  expect_error(
    se_boot(a = c(1, -1), b = c(0, 0), theta = 0, estimator = "MLE"),
    "^a must"
  )
  expect_error(
    design_cumulants(a = c(1, -1), b = c(0, 0), theta = 0), "^a must"
  )
  expect_error(
    delta2(a = c(1, -1), b = c(0, 0), theta = 0, estimator = "MLE"),
    "^a must"
  )
  expect_error(
    delta3(a = c(1, -1), b = c(0, 0), theta = 0, estimator = "MLE"),
    "^a must"
  )
  expect_error(
    risk_diff(
      a = c(1, -1), b = c(0, 0), theta = 0, estimator = "MLE", margin = 1
    ),
    "^a must"
  )
  expect_error(
    se_diagnostic(a = c(1, -1), b = c(0, 0), theta = 0, estimator = "MLE"),
    "^a must"
  )
  expect_error(
    se_select(
      a = c(1, -1), b = c(0, 0), theta = 0, estimator = "MLE", margin = 1
    ),
    "^a must"
  )

  expect_error(
    se_boot(a = c(1, 1), b = c(0, 0), theta = 0), "^estimator must"
  )
  expect_error(
    delta2(a = c(1, 1), b = c(0, 0), theta = 0), "^estimator must"
  )
  expect_error(
    delta3(a = c(1, 1), b = c(0, 0), theta = 0), "^estimator must"
  )
  expect_error(
    risk_diff(a = c(1, 1), b = c(0, 0), theta = 0), "^estimator must"
  )
  expect_error(
    se_diagnostic(a = c(1, 1), b = c(0, 0), theta = 0), "^estimator must"
  )
  expect_error(
    se_select(a = c(1, 1), b = c(0, 0), theta = 0, margin = 1),
    "^estimator must"
  )
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
