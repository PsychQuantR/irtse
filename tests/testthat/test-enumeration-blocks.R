# Evaluate `f` with the estimating equation split into blocks of `size` rows.
with_equation_block <- function(size, f) {
  testthat::local_mocked_bindings(equation_block = function(J) as.integer(size),
                                  .package = "irtse")
  f()
}

# The value of `f()`, or the message of the error it stops with.
outcome_of <- function(f) {
  tryCatch(f(), error = function(e) paste("error:", conditionMessage(e)))
}

# Requirement: Limits of exact enumeration
# Scenario: Root solving in blocks of abilities
test_that("The estimating equation in blocks gives the same doubles and the same error", {
  theta <- seq(-2, 2, length.out = 11)
  a <- c(0.7, 1.3, 2.1)
  b <- c(-0.4, 0.2, 1.1)
  totals <- c(0.5, 1.2, 2.9, 0.1, 3.4, 1.7, 2.2, 0.9, 3.9, 1.1, 2.6)
  calls <- list(
    mle = function() irtse:::estimate_equation(theta, a, b, "MLE", totals),
    wle = function() irtse:::estimate_equation(theta, a, b, "WLE", totals),
    scalar_total = function() irtse:::estimate_equation(theta, a, b, "WLE", 0),
    no_rows = function() irtse:::estimate_equation(numeric(0), a, b, "MLE", numeric(0)),
    # The first row fails only the check of the derivative and the second the
    # check of the value. Evaluated at once the value check comes first, so the
    # error names the estimating equation; a block of one row must not report
    # the first row's derivative instead.
    error_order = function() {
      irtse:::estimate_equation(c(1e-308, 1e-308), c(1.5e308, 1.5e308), c(0, 0), "MLE",
                                c(1.7e308, 0))
    },
    # The check of the logits runs inside each block. The first row passes it and
    # fails only the value check; the logit of the second row is not finite. At once
    # the logits are checked first, so the error names a logit; a block of one row
    # reaches the second row's logit before any deferred check runs and must name
    # a logit as well.
    logit_after_value = function() {
      irtse:::estimate_equation(c(1e-308, 2), c(1.5e308, 1.5e308), c(0, 0), "MLE", c(0, 0))
    }
  )
  expect_match(outcome_of(calls$error_order), "estimating equation is not representable",
               fixed = TRUE)
  expect_match(outcome_of(calls$logit_after_value), "a logit is not representable",
               fixed = TRUE)
  for (name in names(calls)) {
    reference <- outcome_of(calls[[name]])
    for (size in c(1, 2, 7)) {
      expect_identical(with_equation_block(size, function() outcome_of(calls[[name]])),
                       reference, info = paste(name, size))
    }
  }
})

# Requirement: Limits of exact enumeration
# Scenario: Root solving in blocks of abilities
test_that("se_boot() and risk_diff() return the same doubles when root solving runs in blocks", {
  a <- rep(c(0.7, 1.4), 4)
  b <- seq(-1.2, 1.5, length.out = 8)
  calls <- function() {
    list(se_boot(a, b, 0.3, "MLE"), se_boot(a, b, 0.3, "WLE"),
         se_boot(c(0.7, 1.3, 2.1), c(-0.4, 0.2, 1.1), 0.35, "WLE"),
         risk_diff(a, b, -0.2, "MLE"), risk_diff(a, b, -0.2, "WLE", margin = 0.5),
         se_boot(1, 0, 0, "WLE"), outcome_of(function() se_boot(1, 0, 0, "MLE")))
  }
  reference <- calls()
  for (size in c(1, 2, 7)) {
    expect_identical(with_equation_block(size, calls), reference, info = size)
  }
})

# Replace the enumeration of attainable totals by a function that fails, so that a
# test sees whether the lower bound refused the design before it.
local_failing_enumeration <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    reachable_points = function(weights) stop("the attainable totals were enumerated"),
    .package = "irtse", .env = env)
}

# Requirement: Limits of exact enumeration
# Scenario: Designs refused by a lower bound of the atom count
# spec-example: irtse-standard-errors / atom count lower bound
test_that("A lower bound of the atom count above the limit refuses the design before enumeration", {
  local_failing_enumeration()
  rows <- list(
    list(call = function() se_boot(rep(1, 100002), rep(0, 100002), 0, "MLE"),
         bound = "at least 100001 atoms", limit = "se_boot() limit of 100000"),
    list(call = function() se_boot(rep(1, 100000), rep(0, 100000), 0, "WLE"),
         bound = "at least 100001 atoms", limit = "se_boot() limit of 100000"),
    list(call = function() risk_diff(rep(1, 2002), rep(0, 2002), 0, "MLE"),
         bound = "at least 2001 atoms", limit = "risk_diff() limit of 2000"),
    list(call = function() risk_diff(rep(1, 2000), rep(0, 2000), 0, "WLE"),
         bound = "at least 2001 atoms", limit = "risk_diff() limit of 2000")
  )
  for (row in rows) {
    message <- outcome_of(row$call)
    expect_match(message, "exact enumeration is not available", fixed = TRUE)
    expect_match(message, row$bound, fixed = TRUE)
    expect_match(message, row$limit, fixed = TRUE)
  }
})

# Requirement: Limits of exact enumeration
# Scenario: Designs refused by a lower bound of the atom count
test_that("A design whose lower bound equals the limit is enumerated", {
  local_failing_enumeration()
  # The MLE removes the two extreme totals, so its bound is J - 1; the WLE keeps
  # them, so its bound is J + 1. At these sizes the bound equals the limit, so
  # the design reaches the enumeration (here replaced by a function that fails).
  calls <- list(
    function() se_boot(rep(1, 100001), rep(0, 100001), 0, "MLE"),
    function() se_boot(rep(1, 99999), rep(0, 99999), 0, "WLE"),
    function() risk_diff(rep(1, 2001), rep(0, 2001), 0, "MLE"),
    function() risk_diff(rep(1, 1999), rep(0, 1999), 0, "WLE")
  )
  for (call in calls) {
    expect_identical(outcome_of(call), "error: the attainable totals were enumerated")
  }
})
