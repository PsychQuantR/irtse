# Requirement: Design cumulants
test_that("Cumulants of rational designs", {
  cases <- list(
    list(
      a = c(1, 1, 1), b = c(-log(4), -log(4), log(4)), theta = 0,
      expected = c(k2 = 4/25, k3 = -4/125, k4 = 4/625, k5 = 92/3125)
    ),
    list(
      a = c(1, 1, 3), b = c(-log(7/6), -log(7/6), -log(12)/3),
      theta = 0,
      expected = c(
        k2 = 64/169, k3 = -1216/2197, k4 = 29104/28561, k5 = -19840/28561
      )
    ),
    list(
      a = c(1, 1), b = c(0, log(4)), theta = 0,
      expected = c(
        k2 = 41/200, k3 = 6/125, k4 = -593/10000, k5 = -138/3125
      )
    ),
    list(
      a = 1, b = 0, theta = 0,
      expected = c(k2 = 1/4, k3 = 0, k4 = -1/8, k5 = 0, k6 = 1/4, k7 = 0)
    ),
    list(
      a = c(1, 1), b = c(log(5), -log(5)), theta = 0,
      expected = c(k2 = 5/36, k4 = 5/216, k6 = -115/972)
    )
  )

  for (case in cases) {
    cumulants <- design_cumulants(a = case$a, b = case$b, theta = case$theta)
    expect_type(cumulants, "double")
    expect_named(cumulants, c("k2", "k3", "k4", "k5", "k6", "k7"))

    for (component in names(case$expected)) {
      expected <- case$expected[[component]]
      if (expected == 0) {
        expect_lt(abs(cumulants[[component]] - expected), 1e-14)
      } else {
        expect_equal(cumulants[[component]], expected, tolerance = 1e-10)
      }
    }
  }
})

# Requirement: Second-order coefficient
test_that("Exact values of Delta_2", {
  cases <- list(
    list(a = 1, b = 0, theta = qlogis(3/10), estimator = "MLE",
         expected = -48800/3087),
    list(a = 1, b = 0, theta = qlogis(1/5), estimator = "MLE",
         expected = -43425/512),
    list(a = 1, b = 0, theta = qlogis(2/5), estimator = "MLE",
         expected = -1475/576),
    list(a = 1, b = 0, theta = 0, estimator = "MLE", expected = 0),
    list(a = c(1, 1, 1), b = c(-log(4), -log(4), log(4)), theta = 0,
         estimator = "MLE", expected = 2175/512),
    list(a = c(1, 1, 3), b = c(-log(7/6), -log(7/6), -log(12)/3),
         theta = 0, estimator = "MLE", expected = 19410495/4194304),
    list(a = c(1, 1), b = c(0, log(4)), theta = 0, estimator = "MLE",
         expected = -1071573811200/194754273881),
    list(a = c(1, 1), b = c(log(5), -log(5)), theta = 0,
         estimator = "MLE", expected = 0),
    list(a = c(1, 1, 1), b = c(log(7), log(3/5), log(3/5)), theta = 0,
         estimator = "MLE", expected = 80876309760/94931877133),
    list(a = c(1, 1, 1), b = c(log(7), log(3/5), log(3/5)), theta = 0,
         estimator = "WLE", expected = -8286245568/94931877133),
    list(a = 1, b = 0, theta = qlogis(1/5), estimator = "WLE",
         expected = -15075/512),
    list(a = 1, b = 0, theta = 0, estimator = "WLE", expected = 0)
  )

  for (case in cases) {
    coefficient <- delta2(
      a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
    )
    if (case$expected == 0) {
      expect_lt(abs(coefficient - case$expected), 1e-12)
    } else {
      expect_equal(coefficient, case$expected, tolerance = 1e-10)
    }
  }
})

# Requirement: Second-order coefficient
test_that("Coefficients do not depend on the number of copies of a block", {
  a <- c(1, 1, 1)
  b <- c(-log(4), -log(4), log(4))

  for (estimator in c("MLE", "WLE")) {
    expect_equal(
      delta2(a = a, b = b, theta = 0, estimator = estimator),
      delta2(a = rep(a, 4), b = rep(b, 4), theta = 0, estimator = estimator),
      tolerance = 1e-12
    )
  }
})

# Requirement: Second-order coefficient
test_that("Scaling under a common discrimination", {
  expect_equal(
    delta2(
      a = rep(2, 3), b = c(-0.4, 0.1, 0.9), theta = 0.3, estimator = "MLE"
    ),
    delta2(
      a = rep(1, 3), b = 2 * c(-0.4, 0.1, 0.9), theta = 0.6,
      estimator = "MLE"
    ) / 4,
    tolerance = 1e-10
  )
})

# Requirement: Third-order coefficient
test_that("Exact values of Delta_3", {
  cases <- list(
    list(a = 1, b = 0, theta = 0, expected = -121/3),
    list(a = c(1, 1), b = c(log(5), -log(5)), theta = 0,
         expected = 110889/625),
    list(a = 1, b = 0, theta = qlogis(3/10), expected = -369968225/583443),
    list(a = 1, b = 0, theta = qlogis(1/5), expected = -914779225/196608)
  )

  for (case in cases) {
    expect_equal(
      delta3(a = case$a, b = case$b, theta = case$theta, estimator = "MLE"),
      case$expected,
      tolerance = 1e-10
    )
  }
})

# Requirement: Third-order coefficient
test_that("Scaling of Delta_3 under a common discrimination", {
  expect_equal(
    delta3(
      a = rep(2, 3), b = c(-0.4, 0.1, 0.9), theta = 0.3, estimator = "MLE"
    ),
    delta3(
      a = rep(1, 3), b = 2 * c(-0.4, 0.1, 0.9), theta = 0.6,
      estimator = "MLE"
    ) / 4,
    tolerance = 1e-10
  )
})

# Requirement: Third-order coefficient
test_that("No third-order coefficient for the WLE", {
  expect_message(
    coefficient <- delta3(a = 1, b = 0, theta = 0, estimator = "WLE"),
    "WLE"
  )
  expect_true(is.na(coefficient))
  expect_type(coefficient, "double")
})
