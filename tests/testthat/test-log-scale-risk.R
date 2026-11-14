risk_error_message <- function(expr) {
  tryCatch({ force(expr); "" }, error = conditionMessage)
}

upper_risk_outcome_conforms <- function(expr) {
  value <- tryCatch(force(expr), error = identity)
  if (inherits(value, "error")) {
    return(grepl("scaled risk smse_info is not representable",
                 conditionMessage(value), fixed = TRUE))
  }
  is.numeric(value) && length(value) == 1L && is.finite(value) &&
    abs(value / .Machine$double.xmax - 1) <= 1e-12
}

# Requirement: Numerical range
# Scenario: A positive risk below the smallest positive double
test_that("Positive risks below the smallest double name smse_boot", {
  for (d in c(1e-9, 3e-9)) {
    expect_error(risk_diff(rep(1000, 3), c(-1.44, 0, 1.44), -0.72 + d, "MLE"),
                 "scaled risk smse_boot is not representable", fixed = TRUE)
  }
  expect_error(risk_diff(rep(512, 3), c(-2.8125, 0, 2.8125),
                         -1.40625 + 2^-30, "MLE"),
               "scaled risk smse_boot is not representable", fixed = TRUE)
})

# Requirement: Numerical range
# Scenario: A risk representable only as a subnormal number
test_that("Subnormal risks are accurate to one smallest-double spacing", {
  displacements <- c(5e-9, 1e-8, 1e-7)
  units <- c(1.599, 6.397, 639.73)
  for (i in seq_along(displacements)) {
    value <- risk_diff(rep(1000, 3), c(-1.44, 0, 1.44),
                       -0.72 + displacements[[i]], "MLE")$smse_boot
    expect_lte(abs(value / 2^-1074 - units[[i]]), 1)
  }
})

# Requirement: Numerical range
# Scenario: A risk representable only as a subnormal number
test_that("Upper subnormal risks meet the relative accuracy bound", {
  for (K in c(712, 714, 716)) {
    r <- K / 512
    got <- risk_diff(rep(512, 3), c(-2 * r, 0, 2 * r), 2 * r, "MLE")$smse_boot
    u <- r * exp(-K / 2)
    # The operands stay normal until this final multiplication.
    want <- (2 * u) * (6 * u)
    expect_lte(abs(got / want - 1), 1e-12)
  }
})

# Requirement: Numerical range
# Scenario: A risk representable only as a subnormal number
test_that("An exact subnormal power meets the relative accuracy bound", {
  got <- irtse:::scaled_risk_sum(0, 2^-512, what = "smse_boot")
  expect_lte(abs(got / 2^-1024 - 1), 1e-12)
})

# Requirement: Numerical range
# Scenario: A risk at the upper end of the double range
test_that("Risks outside the upper boundary band have one permitted outcome", {
  B0 <- 1042 * log(2) / 512
  B <- B0 * (1 - 1e-9)
  got <- risk_diff(rep(512, 2), c(-B, B), 0, "MLE")$smse_info
  want <- exp(512 * B - 2 * log(512))
  expect_lte(abs(got / want - 1), 1e-12)
  B <- B0 * (1 + 1e-9)
  expect_error(risk_diff(rep(512, 2), c(-B, B), 0, "MLE"),
               "scaled risk smse_info is not representable", fixed = TRUE)
})

# Requirement: Numerical range
# Scenario: A risk at the upper end of the double range
test_that("Risks inside the upper boundary band have two permitted outcomes", {
  B0 <- 1042 * log(2) / 512
  expect_true(upper_risk_outcome_conforms(
    risk_diff(rep(512, 2), c(-B0, B0), 0, "MLE")$smse_info
  ))
  expect_true(upper_risk_outcome_conforms(
    irtse:::scaled_risk_sum(0, 2^512, what = "smse_info")
  ))
})

# Requirement: Numerical range
# Scenario: A risk whose terms leave the double range before they are summed
test_that("A squared report may overflow while its weighted risk is finite", {
  expect_lt(abs(risk_diff(rep(0.7092, 3), c(-1000, 1000, 1000),
                          1000, "MLE")$smse_info /
                  (exp(709.2) / (2 * sqrt(2) * 0.7092^2)) - 1), 1e-8)
})

# Requirement: Numerical range
# Scenario: A risk whose terms leave the double range before they are summed
test_that("An overflowing complete risk names smse_info", {
  expect_error(risk_diff(rep(0.711, 3), c(-1000, 1000, 1000), 1000, "MLE"),
               "scaled risk smse_info is not representable", fixed = TRUE)
})

# Requirement: Numerical range
# Scenario: Extreme abilities and discriminations
test_that("Unscaled zero score and information do not stop the root search", {
  expect_lt(abs(se_boot(rep(1000, 3), c(-1, 1, 1), 0, "MLE") /
                  ((1 + log(2) / 2000) * sqrt(2) * exp(-500)) - 1), 1e-6)
})

# Requirement: Numerical range
# Scenario: Terms of a risk whose weight or partial product is outside the double range
test_that("Risk summation retains a term with an underflowing weight", {
  expect_lte(abs(irtse:::scaled_risk_sum(c(0, -750), c(0, 1000),
                                       what = "smse_boot") / 2^-1074 - 3849.05), 1)
})

# Requirement: Numerical range
# Scenario: Terms of a risk whose weight or partial product is outside the double range
test_that("Risk summation retains a term with an underflowing partial product", {
  expect_lte(abs(irtse:::scaled_risk_sum(-740, 1e-100, 1e100,
                                       what = "diff") / 2^-1074 - 84.78), 1)
})

# Requirement: Numerical range
# Scenario: A risk difference whose terms cancel
test_that("Exactly cancelling risk differences return zero on both scales", {
  for (log_weight in c(0, -800)) {
    # Catching here lets both expectations run when the helper is absent.
    value <- tryCatch(irtse:::scaled_risk_sum(rep(log_weight, 2), c(1, 1),
                                             c(1, -1), what = "diff"),
                      error = conditionMessage)
    expect_identical(value, 0)
  }
  expect_error(irtse:::scaled_risk_sum(c(-800, -800), c(1, 1), c(2, -1),
                                      what = "diff"),
               "scaled risk diff is not representable", fixed = TRUE)
})

# Requirement: Numerical range
# Scenario: Inputs of the risk summation that are not finite
test_that("Risk summation rejects every nonfinite input before dropping zeros", {
  for (argument in c("logs", "first", "second")) {
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      args <- list(logs = c(0, 0), first = c(1, 1), second = c(1, 1),
                   what = "smse_info")
      args[[argument]][[2]] <- value
      expect_match(risk_error_message(do.call(irtse:::scaled_risk_sum, args)),
                   "scaled risk smse_info is not representable", fixed = TRUE)
    }
  }
  expect_error(irtse:::scaled_risk_sum(c(0, -Inf), c(0, 0), what = "smse_info"),
               "scaled risk smse_info is not representable", fixed = TRUE)
})

# Requirement: Numerical range
# Scenarios: A positive risk below the smallest positive double;
# A risk whose terms leave the double range before they are summed
test_that("Each risk field is named for underflow and overflow", {
  for (field in c("smse_info", "smse_boot", "diff", "smse_select")) {
    for (log_weight in c(-745, 710)) {
      expect_match(risk_error_message(irtse:::scaled_risk_sum(log_weight, 1, what = field)),
                   paste("scaled risk", field,
                         "is not representable in double precision."), fixed = TRUE)
    }
  }
})

# Requirement: Numerical range
# Scenario: A risk whose reports agree with the target
test_that("Zero factors give exact zero without exponentiating the weights", {
  value <- tryCatch(irtse:::scaled_risk_sum(c(-800, 710), c(0, 1), c(1, 0),
                                           what = "diff"),
                    error = conditionMessage)
  expect_identical(value, 0)
})

# Requirement: Numerical range
# Scenario: The four fields of risk_diff() use the risk summation
test_that("All four risk fields use the shared summation and specified factors", {
  skip_if_not_installed("testthat", "3.1.7")
  args <- list(a = rep(c(1, 2), 3), b = rep(c(-1/2, 1), 3), theta = 0.3,
               estimator = "MLE", margin = 0)
  ordinary <- do.call(risk_diff, args)
  recorded <- list()
  sentinels <- c(smse_info = 11, smse_boot = 22, diff = 33, smse_select = 44)
  testthat::local_mocked_bindings(scaled_risk_sum = function(logs, first,
                                                           second = first, what) {
    recorded[[length(recorded) + 1L]] <<- list(what = what, logs = logs,
                                             first = first, second = second)
    unname(sentinels[[what]])
  }, .package = "irtse")
  result <- do.call(risk_diff, args)
  for (field in names(sentinels)) {
    expect_identical(result[[field]], unname(sentinels[[field]]))
  }
  expect_length(recorded, 4L)
  names(recorded) <- vapply(recorded, `[[`, "", "what")
  for (field in c("smse_info", "smse_boot", "smse_select")) {
    expect_identical(recorded[[field]]$first, recorded[[field]]$second)
    expect_identical(recorded[[field]]$logs, recorded$diff$logs)
  }
  A <- recorded$smse_info$first + ordinary$h_true
  B <- recorded$smse_boot$first + ordinary$h_true
  expect_equal(recorded$diff$first, A - B, tolerance = 1e-12)
  expect_equal(recorded$diff$second, A + B - 2 * ordinary$h_true, tolerance = 1e-12)
  product <- with(recorded$diff, sum(exp(logs) * first * second))
  # Relational checks against the unmocked call are required by this Scenario.
  expect_lt(abs(product / ordinary$diff - 1), 1e-10)
  expect_lt(abs(product / (ordinary$smse_info - ordinary$smse_boot) - 1), 1e-10)
})
