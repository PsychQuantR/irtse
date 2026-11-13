# Requirement: Exact scaled risks
test_that("Six identical Rasch items at their difficulty", {
  risk <- risk_diff(
    a = rep(1, 6), b = rep(0, 6), theta = 0, estimator = "MLE"
  )
  h0 <- sqrt((36 * log(5)^2 + 90 * log(2)^2) / 31)

  expect_type(risk, "list")
  expect_named(
    risk, c("J", "h_true", "smse_info", "smse_boot", "diff"),
    ignore.order = TRUE
  )
  expect_true(risk$J == 6)
  expect_lt(abs(risk$smse_info - 6960287991367/100000000000000), 1e-12)
  expect_true(round(risk$smse_info, 4) == 696/10000)
  expect_lt(abs(risk$smse_boot - 9724843390404/100000000000000), 1e-12)
  expect_true(round(risk$smse_boot, 4) == 972/10000)
  expect_equal(risk$h_true, h0, tolerance = 1e-12)
  expect_lt(abs(risk$diff - (risk$smse_info - risk$smse_boot)), 1e-14)
})

# Requirement: Exact scaled risks
test_that("The plug-in exact-risk rule can be worse than both fixed reports", {
  a <- rep(1, 6)
  b <- rep(0, 6)
  truth <- risk_diff(a = a, b = b, theta = 0, estimator = "MLE")
  totals <- 1:5
  estimates <- log(totals / (6 - totals))
  weights <- choose(6, totals) / 62

  choose_info <- vapply(estimates, function(estimate) {
    risk_diff(a = a, b = b, theta = estimate, estimator = "MLE")$diff < 0
  }, logical(1))
  info_reports <- vapply(estimates, function(estimate) {
    sqrt(6) * se_info(a = a, b = b, theta = estimate)
  }, numeric(1))
  boot_reports <- vapply(estimates, function(estimate) {
    sqrt(6) * se_boot(a = a, b = b, theta = estimate, estimator = "MLE")
  }, numeric(1))
  reports <- ifelse(choose_info, info_reports, boot_reports)
  rule_risk <- sum(weights * (reports - truth$h_true)^2)

  expect_true(round(rule_risk, 4) == 1004/10000)
  expect_lt(abs(rule_risk - 10036670224707/100000000000000), 1e-12)
  expect_gt(rule_risk, truth$smse_info)
  expect_gt(rule_risk, truth$smse_boot)
  expect_identical(totals[choose_info], 3L)
})

# Requirement: Exact scaled risks
test_that("Sign of the risk difference along a repeated Rasch block", {
  cases <- data.frame(
    J = c(3, 6, 9, 12, 15, 18, 84, 87),
    expected_sign = c(1, -1, 1, 1, 1, -1, -1, 1),
    expected_diff = c(
      NA, NA, NA, NA, NA, NA,
      -145321848/100000000000000, 83070903/100000000000000
    )
  )

  for (i in seq_len(nrow(cases))) {
    J <- cases$J[[i]]
    m <- J / 3
    risk <- risk_diff(
      a = rep(1, J), b = c(rep(0, 2 * m), rep(2, m)), theta = 0.75,
      estimator = "MLE"
    )
    expect_true(sign(risk$diff) == cases$expected_sign[[i]])
    if (!is.na(cases$expected_diff[[i]])) {
      expect_lt(abs(risk$diff - cases$expected_diff[[i]]), 1e-12)
    }
  }
})

# Requirement: Exact scaled risks
test_that("WLE risks on a repeated Rasch block", {
  risk <- risk_diff(
    a = rep(1, 8), b = rep(c(0, 2), 4), theta = 0.5, estimator = "WLE"
  )

  expect_equal(
    risk$smse_info, 798727481557971/10000000000000000, tolerance = 1e-10
  )
  expect_equal(
    risk$smse_boot, 114705194375896/10000000000000000, tolerance = 1e-10
  )
})

# Requirement: Exact scaled risks
test_that("WLE risk difference for one hundred identical items", {
  risk <- risk_diff(
    a = rep(1, 100), b = rep(0, 100), theta = log(1/4), estimator = "WLE"
  )

  expect_equal(100^2 * risk$diff, -456071723352/10000000000, tolerance = 1e-9)
})
