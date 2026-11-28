# Requirement: Estimators on the sufficient statistic
# spec-example: irtse-standard-errors / estimates and root checks
test_that("MLE roots for identical Rasch items", {
  h0 <- sqrt((36 * log(5)^2 + 90 * log(2)^2) / 31)

  expect_equal(
    se_boot(a = rep(1, 6), b = rep(0, 6), theta = 0, estimator = "MLE"),
    h0 / sqrt(6),
    tolerance = 1e-12
  )
})

# Requirement: Estimators on the sufficient statistic
# spec-example: irtse-standard-errors / estimates and root checks
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
# spec-example: irtse-standard-errors / exact laws and identified designs
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
# spec-example: irtse-standard-errors / exact laws and identified designs
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
# spec-example: irtse-standard-errors / exact laws and identified designs
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
# spec-example: irtse-standard-errors / information se for six rasch items
test_that("Six identical Rasch items at their difficulty", {
  expect_equal(
    se_info(a = rep(1, 6), b = rep(0, 6), theta = 0),
    sqrt(2/3),
    tolerance = 1e-12
  )
})

# Requirement: Exact parametric-bootstrap standard error
# spec-example: irtse-standard-errors / bootstrap se and scaled target
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
# spec-example: irtse-standard-errors / lattice enumeration limits
test_that("Lattice limits", {
  # An unmatched error must not prevent the remaining limit checks from running.
  check_error_message <- function(code) {
    tryCatch(force(code), simpleError = function(error) {
      fail(conditionMessage(error))
    })
  }

  check_error_message(expect_error(
    se_boot(c(1, 999999/1000), c(0, 0), 0, "MLE"),
    "lattice width", fixed = TRUE
  ))
  check_error_message(expect_error(
    se_boot(c(1, 999999/1000), c(0, 0), 0, "MLE"),
    "1000999", fixed = TRUE
  ))
  check_error_message(expect_error(
    se_boot(c(1, 999999/1000), c(0, 0), 0, "MLE"),
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
# spec-example: irtse-standard-errors / invalid inputs and the argument named in the error
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
# spec-example: irtse-standard-errors / invalid inputs and the argument named in the error
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
# spec-example: irtse-standard-errors / readme statements and range
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
  expect_no_match(readme, "working name", fixed = TRUE)
  expect_match(readme, "catR::semTheta", fixed = TRUE)
})

# Requirement: Estimators on the sufficient statistic
# spec-example: irtse-standard-errors / estimates and root checks
test_that("Saturated probabilities do not stop the MLE at a wrong root", {
  root <- irtse:::estimate_atoms(1, c(1, 1), c(-40, 42), "MLE")
  expect_lt(abs(root - 1), 1e-12)
  # One atom (T = 1, estimate 1): h_true is 0 and the information report at 1
  # is sqrt(2 / I(1)), I(1) = 2 plogis(41) plogis(-41). Closed form written
  # independently of the package; the wrong root 0 would give 4.15e17.
  risk <- risk_diff(c(1, 1), c(-40, 42), 1, "MLE")$smse_info
  expected <- exp(41) * (1 + exp(-41))^2
  expect_true(is.finite(risk))
  expect_lt(abs(risk / expected - 1), 1e-10)
})

# Requirement: Estimators on the sufficient statistic
# spec-example: irtse-standard-errors / estimates and root checks
test_that("Shifted WLE designs", {
  shifted <- tryCatch(se_boot(1, 1000, 1000, "WLE"), error = identity)
  expect_type(shifted, "double")
  if (is.numeric(shifted)) expect_lt(abs(shifted / log(3) - 1), 1e-12)
  steep <- tryCatch(se_boot(rep(1000, 2), c(0, 0), 0, "WLE"), error = identity)
  expect_true(is.numeric(steep) && is.finite(steep) && steep > 0)
})

# Requirement: Estimators on the sufficient statistic
# spec-example: irtse-standard-errors / estimates and root checks
test_that("Two WLE roots inside one grid cell", {
  expect_error(
    se_boot(c(40, 41, 3.480672), c(0, 0.1, 0), 0, "WLE"),
    "T = 0", fixed = TRUE
  )
})

# Requirement: Numerical range
# spec-example: irtse-standard-errors / numerical range
test_that("Extreme abilities and discriminations", {
  cases <- list(
    list(call = function() se_boot(rep(1, 6), rep(0, 6), 400, "MLE"),
         expected = sqrt(5/2) * log(5/2) * exp(-200), tolerance = 1e-6,
         must_return = TRUE),
    list(call = function() se_info(rep(1, 6), rep(0, 6), 600),
         expected = exp(300) / sqrt(6), tolerance = 1e-10, must_return = TRUE),
    list(call = function() delta2(1, 0, 300, "MLE"), negative = TRUE),
    list(call = function() delta3(1, 0, 75, "MLE"), negative = TRUE,
         must_return = TRUE),
    list(call = function() delta2(1e-3, 0, 0, "MLE"), expected = 0,
         must_return = TRUE),
    list(call = function() delta3(1e-3, 0, 0, "MLE"),
         expected = -121 / (3 * 1e-6), tolerance = 1e-10, must_return = TRUE),
    list(call = function() se_select(rep(1, 4), rep(0, 4), 300, "MLE", 0),
         is_list = TRUE),
    list(call = function() se_diagnostic(rep(1, 6), rep(0, 6), 300,
                                        "WLE")$sign_delta2, is_sign = TRUE)
  )
  for (i in seq_along(cases)) {
    case <- cases[[i]]
    if (i %in% c(3L, 7L, 8L)) {
      expect_error(case$call(), "representable", fixed = TRUE)
      next
    }
    value <- tryCatch(case$call(), error = identity)
    if (inherits(value, "error")) {
      # Only the rows whose spec entry offers the alternative may stop.
      expect_false(isTRUE(case$must_return),
                   info = paste("Numerical range row", i, "must return a value:",
                                conditionMessage(value)))
      expect_match(conditionMessage(value), "representable", fixed = TRUE,
                   info = paste("Numerical range row", i))
    } else if (isTRUE(case$is_list)) {
      expect_type(value, "list")
      reported <- unlist(value[vapply(value, is.numeric, logical(1))])
      expect_false(any(is.nan(reported) | is.infinite(reported)))
    } else {
      expect_true(length(value) == 1L && is.finite(value),
                  info = paste("Numerical range row", i))
      if (i %in% c(1L, 2L)) expect_gt(value, 0)
      if (isTRUE(case$negative)) expect_lt(value, 0)
      if (isTRUE(case$is_sign)) expect_true(value %in% c(-1, 0, 1))
      if (!is.null(case$expected)) {
        if (case$expected == 0) expect_identical(value, 0) else {
          expect_lt(abs(value / case$expected - 1), case$tolerance)
        }
      }
    }
  }
})

# Requirement: Estimators on the sufficient statistic
test_that("Per-T roots for identical Rasch items", {
  for (estimator in c("MLE", "WLE")) {
    totals <- if (estimator == "MLE") 1:5 else 0:6
    expected <- if (estimator == "MLE") log(totals / (6 - totals)) else {
      qlogis((totals + 1/2) / 7)
    }
    roots <- irtse:::estimate_atoms(totals, rep(1, 6), rep(0, 6), estimator)
    for (i in seq_along(totals)) expect_lt(abs(roots[[i]] - expected[[i]]), 1e-12)
  }
})

# Requirement: Estimators on the sufficient statistic
test_that("Per-T asymmetric roots satisfy an independent stable residual", {
  # Increment/decrement the IEEE-754 encoding to obtain neighbouring doubles.
  neighbour <- function(x, up) {
    if (x == 0) return(if (up) 2^-1074 else -2^-1074)
    bytes <- as.integer(writeBin(x, raw(), size = 8, endian = "little"))
    carry <- if (up == (x > 0)) 1L else -1L
    for (j in 1:8) {
      next_byte <- bytes[[j]] + carry
      bytes[[j]] <- next_byte %% 256L
      if (next_byte >= 0L && next_byte <= 255L) break
    }
    readBin(as.raw(bytes), "double", n = 1, size = 8, endian = "little")
  }
  a <- c(1, 2)
  b <- c(-0.5, 0.7)
  residual <- function(t, total, estimator) {
    eta <- a * (t - b)
    positive <- eta > 0
    score <- (total - sum(a[positive])) +
      sum(a[positive] * plogis(-eta[positive])) -
      sum(a[!positive] * plogis(eta[!positive]))
    if (estimator == "MLE") return(score)
    log_w <- 2 * log(a) + plogis(eta, log.p = TRUE) +
      plogis(-eta, log.p = TRUE)
    w <- exp(log_w - max(log_w))
    score + sum(w * a * (plogis(-eta) - plogis(eta))) / (2 * sum(w))
  }
  for (estimator in c("MLE", "WLE")) {
    totals <- if (estimator == "MLE") 1:2 else 0:3
    roots <- irtse:::estimate_atoms(totals, a, b, estimator)
    for (i in seq_along(totals)) {
      values <- vapply(c(neighbour(roots[[i]], FALSE), roots[[i]],
                         neighbour(roots[[i]], TRUE)), residual, numeric(1),
                       total = totals[[i]], estimator = estimator)
      expect_true(values[[2]] == 0 || sign(values[[1]]) * sign(values[[3]]) <= 0,
                  info = paste(estimator, "T =", totals[[i]]))
    }
  }
})

# Requirement: Estimators on the sufficient statistic
test_that("WLE multiple-root message names T", {
  expect_error(se_boot(c(1, 30), c(0, 0), 0, "WLE"), "T = ", fixed = TRUE)
})

# Requirement: Package metadata and documentation
test_that("Package metadata matches the specification", {
  description <- utils::packageDescription("irtse")
  expect_identical(description$Version, "0.0.0.9000")
  expect_identical(description$License, "MIT + file LICENSE")
  expect_identical(trimws(strsplit(description$Imports, ",")[[1]]), "stats")
})

# Requirement: Package metadata and documentation
# spec-example: irtse-standard-errors / authors and contact
test_that("The package names its two authors and its contact", {
  description <- utils::packageDescription("irtse")
  authors <- eval(parse(text = description$`Authors@R`))
  expect_length(authors, 2L)
  expect_identical(authors[[1]]$given, "Che")
  expect_identical(authors[[1]]$family, "Cheng")
  expect_identical(authors[[1]]$email, "che830621@stat.sinica.edu.tw")
  expect_identical(authors[[1]]$role, c("aut", "cre"))
  expect_identical(authors[[2]]$given, "Hau-Hung")
  expect_identical(authors[[2]]$family, "Yang")
  expect_identical(authors[[2]]$email, "hauhungyang@as.edu.tw")
  expect_identical(authors[[2]]$role, "aut")
  # exactly one person is the maintainer and contact, and it is the first
  expect_identical(which(vapply(authors, function(p) "cre" %in% p$role, TRUE)), 1L)
  holder <- grep("^COPYRIGHT HOLDER:", readLines(system.file("LICENSE", package = "irtse")), value = TRUE)
  expect_identical(holder, "COPYRIGHT HOLDER: Che Cheng and Hau-Hung Yang")
})

# Requirement: Package metadata and documentation
test_that("Package files contain no local paths or private issue references", {
  paths <- shipped_scan_files(system.file(package = "irtse"),
                              testthat::test_path("..", ".."))
  # The selection is not empty: it holds the DESCRIPTION of the package. (A
  # count above zero would not show that here, because in a source tree the
  # files of the R folder are added after the selection.)
  expect_true("DESCRIPTION" %in% basename(paths))
  expect_identical(files_matching(paths, shipped_scan_patterns()), character())
})

# Requirement: Numerical range
test_that("Log information and cumulants survive unrepresentable intermediates", {
  expect_lt(abs(se_info(rep(1, 6), rep(0, 6), 800) /
                  (exp(400) / sqrt(6)) - 1), 1e-10)
})

# Requirement: Numerical range
test_that("A saturated item does not erase representable coefficients", {
  theta <- qlogis(1/5)
  # Only the first of two items contributes: k_r is halved, multiplying
  # Delta_2 by 2^3 and Delta_3 by 2^4 relative to the rational one-item case.
  expect_lt(abs(delta2(c(1, 1e3), c(0, 1), theta, "MLE") /
                  (-43425/64) - 1), 1e-10)
  expect_lt(abs(delta2(c(1, 1e3), c(0, 1), theta, "WLE") /
                  (-15075/64) - 1), 1e-10)
  expect_lt(abs(delta3(c(1, 1e3), c(0, 1), theta, "MLE") /
                  (-914779225/12288) - 1), 1e-10)
})

# Requirement: Numerical range
test_that("Symmetric cancellation preserves the remaining small cumulant", {
  small <- plogis(-40) * plogis(40) * (plogis(40) - plogis(-40)) / 3
  first <- design_cumulants(rep(1, 3), c(-1, 40, 1), 0)
  reordered <- design_cumulants(rep(1, 3), c(-1, 1, 40), 0)
  expect_lt(abs(first[["k3"]] / small - 1), 1e-12)
  expect_lt(abs(reordered[["k3"]] / small - 1), 1e-12)
  d2 <- delta2(rep(1, 3), c(-1, 40, 1), 0, "MLE")
  expect_true(is.finite(d2) && d2 < 0)
  expect_lt(abs(d2 / delta2(rep(1, 3), c(-1, 1, 40), 0, "MLE") - 1), 1e-12)
})

# Requirement: Numerical range; Exact scaled risks
test_that("Extreme risk targets preserve a representable standard deviation", {
  expected <- sqrt(6) * sqrt(5/2) * log(5/2) * exp(-400)
  value <- risk_diff(rep(1, 6), rep(0, 6), 800, "MLE")$h_true
  expect_gt(value, 0)
  expect_lt(abs(value / expected - 1), 1e-10)
})

# Requirement: Numerical range; Estimators on the sufficient statistic
# spec-example: irtse-standard-errors / numerical range
test_that("Underflowed score tails cannot certify a false MLE root", {
  expect_error(se_boot(c(1, 2), c(-1000, 1000), 0, "MLE"),
               "representable", fixed = TRUE)
})

# Requirement: Estimators on the sufficient statistic
test_that("An interval without a root stops with an error naming T", {
  # No valid design reaches this path through the exported functions (the
  # interval of the requirement always contains a root), so the internal
  # bracketing step is called with an interval that excludes the root.
  expect_error(irtse:::root_brackets(1, c(1, 1), c(0, 0), "MLE", c(5, 6)),
               "MLE root for T = 1", fixed = TRUE)
  expect_error(irtse:::root_brackets(0, 1, 0, "WLE", c(5, 6)),
               "WLE root for T = 0", fixed = TRUE)
})

# Requirement: Numerical range
# spec-example: irtse-standard-errors / numerical range
test_that("A coefficient outside double range stops the functions that return it", {
  # Delta_2 = -(9/8) exp(900) here: its sign is determined but the value is
  # not representable, and se_diagnostic() and se_select() return the value.
  expect_error(delta2(1, 0, 300, "MLE"), "representable", fixed = TRUE)
  expect_error(se_diagnostic(rep(1, 6), rep(0, 6), 300, "WLE"),
               "representable", fixed = TRUE)
  expect_error(se_select(rep(1, 4), rep(0, 4), 300, "MLE", 0),
               "representable", fixed = TRUE)
})


# Requirement: Numerical range
# spec-example: irtse-standard-errors / numerical range
test_that("Extreme range rows require their specified values and errors", {
  expect_error(se_boot(c(1, 2), c(-1000, 1000), 0, "MLE"),
               "representable", fixed = TRUE)
  expect_lt(abs(delta3(1, 0, 75, "MLE") / (-(9689/768) * exp(300)) - 1),
            1e-10)
})

# Requirement: Numerical range; Information-based standard error
# spec-example: irtse-standard-errors / numerical range
test_that("Information does not require representable higher cumulants", {
  expect_lt(abs(se_info(c(1, 1), c(0, 1000), 0) / 2 - 1), 1e-12)
})

# Requirement: Numerical range; Exact scaled risk of the selection rule
# spec-example: irtse-standard-errors / numerical range
test_that("Risk evaluates both MLE coefficients even outside the margin", {
  expect_error(risk_diff(c(1, 2), c(-150, 150), 0, "MLE", margin = 0),
               "representable", fixed = TRUE)
  expect_error(se_select(c(1, 2), c(-150, 150), (150 - log(2)) / 3,
                         "MLE", margin = 0), "representable", fixed = TRUE)
})

# Requirement: Numerical range
# spec-example: irtse-standard-errors / small logits
test_that("Small logits keep their odd cumulants", {
  expect_lt(abs(design_cumulants(1, 0, 1e-18)[["k3"]] / -1.25e-19 - 1), 1e-10)
  expect_lt(abs(delta2(1, 0, 1e-18, "MLE") / -1.4e-35 - 1), 1e-10)
  expect_lt(abs(delta2(1, 0, 1e-18, "WLE") / -5e-36 - 1), 1e-10)
  expect_lt(abs(delta2(1e-3, 0, 1e-17, "MLE") / -1.4e-33 - 1), 1e-10)
  expect_lt(abs(delta2(1e-3, 0, 1e-17, "WLE") / -5e-34 - 1), 1e-10)
  selected <- se_select(rep(1e-3, 2), rep(0, 2), 1e-14, "WLE", margin = 0)
  expect_identical(selected$choice, "info")
  expect_lt(abs(selected$delta2 / -5e-28 - 1), 1e-10)
  expect_lt(abs(se_diagnostic(rep(1, 6), rep(0, 6), 1e-18,
                             "WLE")$sign_delta2 / -1 - 1), 1e-10)
})

# Requirement: Numerical range
# spec-example: irtse-standard-errors / estimates and root checks
test_that("Large common shifts preserve law-based reports", {
  expect_lt(abs(se_boot(1, 1e3, 1e3, "WLE") / log(3) - 1), 1e-12)
  expect_lt(abs(se_boot(1, -1e3, -1e3, "WLE") / log(3) - 1), 1e-12)
  r <- risk_diff(1, 1e3, 1e3, "WLE", margin = 0)
  expect_lt(abs(r$h_true / log(3) - 1), 1e-10)
  expect_lt(abs(r$smse_info / (4 / sqrt(3) - log(3))^2 - 1), 1e-10)
  expect_lt(abs(r$smse_boot / ((sqrt(3) / 2 - 1) * log(3))^2 - 1), 1e-10)
  expect_lt(abs(r$smse_select / (4 / sqrt(3) - log(3))^2 - 1), 1e-10)
  expect_lt(abs(se_select(1, 1e3, 1e3, "WLE", margin = 0)$se_boot /
                  log(3) - 1), 1e-12)
})

# Requirement: Numerical range
test_that("Cancellation between items removes only exact opposite pairs", {
  expected <- plogis(-40) * plogis(40) * (plogis(40) - plogis(-40)) / 3
  for (b in list(c(-1, 40, 1), c(-1, 1, 40))) {
    expect_lt(abs(design_cumulants(rep(1, 3), b, 0)[["k3"]] / expected - 1),
              1e-12)
  }
  value <- tryCatch(design_cumulants(rep(1, 2),
                                    c(1, -(1 + 2^-50)), 0),
                    error = identity)
  if (inherits(value, "error")) {
    expect_match(conditionMessage(value), "not resolved", fixed = TRUE)
  } else {
    expect_true(is.finite(value[["k3"]]) && value[["k3"]] != 0)
  }
})

# Requirement: Input validation at the boundary
# spec-example: irtse-standard-errors / invalid inputs and the argument named in the error
test_that("One-item MLE errors start with the required prefix", {
  expect_error(se_boot(1, 0, 0, "MLE"), "^MLE requires")
})

# Requirement: Documented scope of the selection rule
test_that("Selection scope requires a compact set with no zero", {
  rd_path <- testthat::test_path("..", "..", "man", "se_select.Rd")
  rd <- if (file.exists(rd_path)) readLines(rd_path, warn = FALSE) else {
    as.character(tools::Rd_db("irtse")[["se_select.Rd"]])
  }
  readme <- testthat::test_path("..", "..", "README.md")
  if (!file.exists(readme)) readme <- system.file("README.md", package = "irtse")
  texts <- list(paste(as.character(rd), collapse = " "))
  if (file.exists(readme)) texts <- c(texts, list(readLines(readme, warn = FALSE)))
  for (text in texts) {
    text <- gsub("\\s+", " ", paste(text, collapse = " "))
    expect_match(text, "compact set", fixed = TRUE)
    expect_match(text, "no zero", fixed = TRUE)
  }
})

# Requirement: Package metadata and documentation
test_that("Suggests declares the required testthat version", {
  suggests <- gsub("\\s+", " ", utils::packageDescription("irtse")$Suggests)
  expect_match(suggests, "testthat (>= 3.0.0)", fixed = TRUE)
})

# Requirement: Package metadata and documentation
test_that("All source files exclude local paths and private tracker references", {
  paths <- source_scan_files(system.file(package = "irtse"),
                             testthat::test_path("..", ".."))
  expect_true("DESCRIPTION" %in% basename(paths))
  expect_identical(files_matching(paths, source_scan_patterns()), character())
})

# Requirement: Package metadata and documentation
# A control of the two scans themselves; the specification has no Scenario of
# this name.
test_that("Both scans report a local path planted in a package tree", {
  base <- tempfile("irtse-scan-")
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  put <- function(root, relative, text) {
    path <- file.path(root, relative)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(text, path)
  }
  local_path <- paste0("/", "Users", "/", "someone", "/", "notes.txt")
  report <- function(package_root, source_root) {
    list(shipped = files_matching(shipped_scan_files(package_root, source_root),
                                  shipped_scan_patterns()),
         source = files_matching(source_scan_files(package_root, source_root),
                                 source_scan_patterns()))
  }

  # A source tree. It stands below a folder named .git, which is not part of
  # the path of any file inside the tree. Its own .git folder and the .git
  # file of a fixture below tests hold a local path, as version control does.
  tree <- file.path(base, ".git", "pkg")
  put(tree, "DESCRIPTION", "Package: control")
  put(tree, "README.md", "A control tree.")
  put(tree, "NEWS.md", "No news.")
  put(tree, file.path("R", "code.R"), "f <- function() 1")
  put(tree, file.path("man", "f.Rd"), "\\name{f}")
  put(tree, file.path("tests", "testthat", "test-f.R"), "f()")
  put(tree, file.path(".git", "logs", "HEAD"), local_path)
  put(tree, file.path("tests", "testthat", "fixture", ".git"),
      paste("gitdir:", local_path))
  expect_true("DESCRIPTION" %in% basename(shipped_scan_files(tree, tree)))
  expect_true("DESCRIPTION" %in% basename(source_scan_files(tree, tree)))
  expect_identical(report(tree, tree),
                   list(shipped = character(), source = character()))
  put(tree, file.path("R", "code.R"),
      c("f <- function() 1", paste("# see", local_path)))
  planted <- file.path(tree, "R", "code.R")
  expect_identical(report(tree, tree), list(shipped = planted, source = planted))

  # An installed package as R CMD check lays it out: the folder two levels
  # above the tests, the check folder, has no R folder. Its serialized files
  # can match a pattern by chance and are not read.
  check <- file.path(base, "check")
  installed <- file.path(check, "pkg")
  put(installed, "DESCRIPTION", "Package: control")
  put(installed, "NEWS.md", "No news.")
  put(installed, file.path("R", "pkg.rdb"), local_path)
  put(installed, file.path("help", "aliases.rds"), local_path)
  expect_true("DESCRIPTION" %in% basename(shipped_scan_files(installed, check)))
  expect_true("DESCRIPTION" %in% basename(source_scan_files(installed, check)))
  expect_identical(report(installed, check),
                   list(shipped = character(), source = character()))
  put(installed, "NEWS.md", paste("See", local_path))
  planted <- file.path(installed, "NEWS.md")
  expect_identical(report(installed, check),
                   list(shipped = planted, source = planted))
})

# Requirement: Exact law of the sufficient statistic
# spec-example: irtse-standard-errors / exact laws and identified designs
test_that("Nearly equal discriminations are not merged", {
  d <- 2^-40
  expected <- 2 * d / (1 + (1 + d)^2)
  expect_lt(abs(se_boot(c(1, 1 + d), c(0, 0), 0, "MLE") / expected - 1), 1e-6)
  expect_null(irtse:::lattice_weights(c(1, 1 + d)))
})

# Requirement: Exact law of the sufficient statistic; Limits of exact enumeration
# spec-example: irtse-standard-errors / exact laws and identified designs
test_that("Representation tolerance preserves decimal and earlier lattices", {
  expect_lt(abs(se_boot(c(0.1, 0.2, 0.3), c(0, 0, 0), 0, "WLE") /
                  (10 * se_boot(c(1, 2, 3), c(0, 0, 0), 0, "WLE")) - 1), 1e-10)
  cases <- list(c(0.1, 0.2, 0.3), c(1, 1000001/1000), rep(c(1, 997), 600),
                c(1, 999), c(1, 2), c(1, 2, 3, 1.5), rep(2, 5))
  for (a in cases) expect_false(is.null(irtse:::lattice_weights(a)))
})

# Requirement: Estimators on the sufficient statistic
# spec-example: irtse-standard-errors / estimates and root checks
test_that("Three WLE roots closer together than the grid are counted", {
  B <- acosh(4) + 1e-8
  expect_error(irtse:::estimate_atoms(2, rep(1, 3), c(-100, -B, B), "WLE"),
               "T = 2", fixed = TRUE)
})

# Requirement: Exact scaled risk of the selection rule; Numerical range
test_that("An underflowed no-choice atom still makes selection risk undefined", {
  r <- risk_diff(rep(1, 2), rep(0, 2), 800, "WLE", margin = 0)
  expect_true(is.na(r$smse_select))
  expect_identical(r$p_select[["none"]], 0)
})

# Requirement: Estimators on the sufficient statistic
test_that("Compensated MLE score preserves the saturated asymmetric tail", {
  score <- irtse:::estimate_equation(0, rep(1, 3), c(-40, 0, 0), "MLE", T = 2)$value
  expect_lt(abs(score / (1 / (1 + exp(40))) - 1), 1e-10)
})

# Requirement: Estimators on the sufficient statistic
test_that("The saturated asymmetric MLE root matches its closed form", {
  # With x = exp(t), the T = 2 equation reduces to x^2 - x - 2 exp(-40) = 0.
  # Rationalising x - 1 avoids cancellation in the positive quadratic root.
  expected <- log1p(4 * exp(-40) / (sqrt(1 + 8 * exp(-40)) + 1))
  root <- irtse:::estimate_atoms(2, rep(1, 3), c(-40, 0, 0), "MLE")
  expect_lt(abs(root - expected), 1e-25)
})

# Requirement: Numerical range; Diagnostic without choosing
test_that("An unresolved average stops, and is a zero in the diagnostic search", {
  # Logits 1 and -(1 + 2^-52) are not an exact pair, but their k3 terms round
  # to the same magnitude, so the computed sum is exactly 0.
  expect_error(design_cumulants(rep(1, 2), c(-1, 1 + 2^-52), 0),
               "not resolved", fixed = TRUE)
  factors <- irtse:::coefficient_factors(rep(1, 2), c(-1, 1 + 2^-52), 0, "MLE")
  expect_identical(factors[[1L]], 0)
})

# Requirement: Numerical range; Exact law of the sufficient statistic
test_that("Law-based reports preserve exactly representable common shifts", {
  a <- c(1, 2, 1.5)
  # These binary fractions and their shifted inputs are exactly representable.
  b <- round(c(-0.5, 0.3, 1.2) * 2^30) / 2^30
  theta <- round(0.4 * 2^30) / 2^30
  shift <- 2^9
  for (estimator in c("MLE", "WLE")) {
    boot <- se_boot(a, b, theta, estimator)
    shifted_boot <- se_boot(a, b + shift, theta + shift, estimator)
    expect_lt(abs(shifted_boot / boot - 1), 1e-12)
    risk <- risk_diff(a, b, theta, estimator)
    shifted_risk <- risk_diff(a, b + shift, theta + shift, estimator)
    for (field in c("J", "h_true", "smse_info", "smse_boot", "diff")) {
      expect_lt(abs(shifted_risk[[field]] / risk[[field]] - 1), 1e-12)
    }
    selected <- se_select(a, b, theta, estimator, margin = 0)
    shifted_selected <- se_select(a, b + shift, theta + shift, estimator,
                                  margin = 0)
    expect_lt(abs(shifted_selected$se_boot / selected$se_boot - 1), 1e-12)
  }
})

# Requirement: Numerical range
# spec-example: irtse-standard-errors / inputs outside the supported range
test_that("Delta_2 rejects a discrimination outside the supported range", {
  expect_error(delta2(c(1, 1e308), c(0, 10), 0, "MLE"),
               "^a must")
})

# Requirement: Numerical range
test_that("Information rejects a discrimination outside the supported range", {
  expect_error(se_info(c(1, 1e308), c(0, 10), 0),
               "^a must")
})

# Requirement: Numerical range
test_that("Cumulants reject a discrimination outside the supported range", {
  expect_error(design_cumulants(c(1, 1e308), c(0, 10), 0),
               "^a must")
})

# Requirement: Numerical range
test_that("Bootstrap rejects a discrimination outside the supported range", {
  expect_error(se_boot(c(1, 1e308), c(0, 10), 0, "WLE"),
               "^a must")
})

# Requirement: Numerical range
# spec-example: irtse-standard-errors / numerical range
# spec-example: irtse-standard-errors / inputs outside the supported range
test_that("A nonzero logit whose half underflows is rejected", {
  expect_error(delta2(2^-1074, 0, 1, "MLE"), "^a must")
  expect_error(delta2(1, 0, 2^-1074, "MLE"), "representable", fixed = TRUE)
})

# Requirement: Numerical range
test_that("A representable small logit retains the third cumulant", {
  expect_lt(abs(design_cumulants(1, 0, 2^-1000)[["k3"]] / -2^-1003 - 1), 1e-10)
})

# Requirement: Numerical range; Exact law of the sufficient statistic
test_that("Discriminations inside the lattice tolerance identify one design", {
  exact <- se_boot(c(1, 1), c(0, 0), 0, "MLE")
  expect_identical(exact, 0)
  expect_identical(se_boot(c(1, 1 + 2^-50), c(0, 0), 0, "MLE"), exact)
})

# Requirement: Input validation at the boundary
test_that("Complex inputs name the invalid argument", {
  expect_error(se_info(1 + 1i, 0, 0), "^a must")
  expect_error(se_info(1, 1i, 0), "^b must")
  expect_error(se_info(1, 0, 1i), "^theta must")
  expect_error(se_select(c(1, 1), c(0, 0), 0, "WLE", margin = 1 + 1i),
               "^margin must")
})

# Requirement: Package metadata and documentation
test_that("Local path patterns catch POSIX and forward-slash drive paths", {
  forbidden <- c(paste0("/", "tmp", "/x.R"),
                 paste0("/", "var", "/", "folders", "/x"),
                 paste0("/", "private", "/x"),
                 paste0("/", "srv", "/", "irtse", "/private.txt"),
                 paste0("C", ":", "/", "work", "/x.R"))
  matches <- function(text) {
    any(vapply(additional_local_path_patterns(), function(pattern) {
      length(grepRaw(pattern, charToRaw(text))) > 0L
    }, logical(1)))
  }
  for (text in forbidden) expect_true(matches(text))
  expect_false(matches("https://example.org/a"))
  expect_false(matches(paste0("https://example.org/", "tmp", "/manual.html")))
})

# Requirement: Estimators on the sufficient statistic
test_that("Score tails representable only as subnormal numbers are retained", {
  s <- 1000
  B <- 0.72
  expected_root <- log(2) / (2 * s)
  expected_se <- (B + expected_root) * sqrt(2) * exp(-s * B / 2)
  for (direction in c(-1, 1)) {
    b <- direction * c(B, -B, -B)
    total <- if (direction == -1) s else 2 * s
    root <- irtse:::estimate_atoms(total, rep(s, 3), b, "MLE")
    expect_lt(abs(root / (direction * expected_root) - 1), 1e-6)
    expect_lt(abs(se_boot(rep(s, 3), b, 0, "MLE") / expected_se - 1), 1e-6)
  }
})

# Requirement: Exact law of the sufficient statistic
# spec-example: irtse-standard-errors / exact laws and identified designs
test_that("Identified design with separated difficulties keeps its symmetric root", {
  a <- c(1, 1 + 2^-50)
  b <- c(50, -50)
  result <- risk_diff(a, b, 0, "MLE")
  expect_identical(result$h_true, 0)
  expect_identical(result$smse_boot, 0)
  expect_lt(abs(result$smse_info / (exp(50) * (1 + exp(-50))^2) - 1), 1e-10)
  expect_identical(se_boot(a, b, 0, "MLE"), 0)
})

# Requirement: Exact law of the sufficient statistic
# spec-example: irtse-standard-errors / exact laws and identified designs
test_that("Identified designs give identical risks and rule reports at every atom", {
  a <- c(1, 1 + 2^-50, 1)
  b <- c(50, -50, 0)
  result <- risk_diff(a, b, 0.3, "MLE")
  expected <- risk_diff(c(1, 1, 1), b, 0.3, "MLE")
  for (field in names(expected)) {
    expect_identical(result[[field]], expected[[field]], info = field)
  }
  result <- risk_diff(a, b, 0.3, "MLE", margin = 0)
  expected <- risk_diff(c(1, 1, 1), b, 0.3, "MLE", margin = 0)
  for (field in names(expected)) {
    expect_identical(result[[field]], expected[[field]], info = field)
  }
})

# Requirement: Exact law of the sufficient statistic; Numerical range
test_that("Selection identifies only its law-based bootstrap report", {
  a <- c(1, 1 + 2^-50, 1)
  b <- c(50, -50, 0)
  selected <- se_select(a, b, 0.3, "MLE", margin = 0)
  expect_identical(selected$se_boot, se_boot(c(1, 1, 1), b, 0.3, "MLE"))
  expect_identical(selected$se_info, se_info(a, b, 0.3))
  expect_identical(selected$delta2, delta2(a, b, 0.3, "MLE"))
  expect_identical(selected$delta3, delta3(a, b, 0.3, "MLE"))
})

# Local definitions from the specification, apart from the package: kappa_r and
# k_r of the requirement "Design cumulants", and the formula of the requirement
# "Third-order coefficient".
spec_logistic <- function(eta) {
  ifelse(eta < 0, exp(eta) / (1 + exp(eta)), 1 / (1 + exp(-eta)))
}
spec_cumulants <- function(a, b, theta) {
  P <- spec_logistic(a * (theta - b))
  u <- P * (1 - P)
  kappa <- list(u, u * (1 - 2 * P), u * (1 - 6 * u),
                u * (1 - 2 * P) * (1 - 12 * u), u * (1 - 30 * u + 120 * u^2),
                u * (1 - 2 * P) * (1 - 60 * u + 360 * u^2))
  stats::setNames(vapply(2:7, function(r) mean(a^r * kappa[[r - 1L]]),
                         numeric(1)), paste0("k", 2:7))
}
spec_third_order <- function(k) {
  k2 <- k[["k2"]]; k3 <- k[["k3"]]; k4 <- k[["k4"]]
  k5 <- k[["k5"]]; k6 <- k[["k6"]]; k7 <- k[["k7"]]
  R3 <- 288 * k2^4 * k3 * k7 + 288 * k2^4 * k4 * k6 + 384 * k2^4 * k5^2 -
    5152 * k2^3 * k3^2 * k6 - 16928 * k2^3 * k3 * k4 * k5 -
    2120 * k2^3 * k4^3 + 49552 * k2^2 * k3^3 * k5 +
    93388 * k2^2 * k3^2 * k4^2 - 279274 * k2 * k3^4 * k4 + 169263 * k3^6
  -R3 / (768 * k2^10)
}

# exp(-730) / 2, the value that the first row of the Example "errors from
# intermediate values" names, as a multiple of 2^-1074 (computed with 80
# significant digits; the nearest double is 933713 multiples).
k3_intermediate_multiple <- 933713.328055

# Whether a value that a call of that Example returned conforms to the Scenario
# "Errors from intermediate values": finite, not 0, and within relative 1e-8 of
# the stated value, or, where the stated value is a subnormal number (given as
# its multiple of 2^-1074), within two multiples of 2^-1074 of it.
intermediate_value_conforms <- function(value, stated, subnormal = FALSE) {
  if (!(length(value) == 1L && is.finite(value) && value != 0)) return(FALSE)
  if (subnormal) abs(value / 2^-1074 - stated) <= 2 else abs(value / stated - 1) < 1e-8
}

# Requirement: Numerical range
# Scenario: Errors from intermediate values
# spec-example: irtse-standard-errors / errors from intermediate values
test_that("Errors from intermediate values may stop an otherwise representable result", {
  k <- spec_cumulants(c(1, 1000), c(730, 0), 0)
  cases <- list(
    list(call = function() design_cumulants(c(1, 1000), c(730, 0), 0)[["k3"]],
         stated = k3_intermediate_multiple, subnormal = TRUE, reason = "cumulant ratio"),
    list(call = function() delta3(c(1, 1000), c(730, 0), 0, "MLE"),
         stated = spec_third_order(k), subnormal = FALSE, reason = "cumulant ratio")
  )
  for (case in cases) {
    value <- tryCatch(case$call(), error = identity)
    if (inherits(value, "error")) {
      expect_match(conditionMessage(value), "representable", fixed = TRUE)
      expect_match(conditionMessage(value), case$reason, fixed = TRUE)
    } else {
      expect_true(intermediate_value_conforms(value, case$stated, case$subnormal))
    }
  }
  # Where the package stops on these calls, the comparison above does not run;
  # it is shown here on planted values. It accepts the double nearest to
  # exp(-730) / 2 and rejects one three multiples of 2^-1074 away and 0; it
  # accepts a Delta_3 a relative 1e-9 off and rejects one 1e-7 off.
  expect_true(intermediate_value_conforms(933713 * 2^-1074, k3_intermediate_multiple, TRUE))
  expect_false(intermediate_value_conforms(933716 * 2^-1074, k3_intermediate_multiple, TRUE))
  expect_false(intermediate_value_conforms(0, k3_intermediate_multiple, TRUE))
  d3 <- spec_third_order(k)
  expect_true(intermediate_value_conforms(d3 * (1 + 1e-9), d3))
  expect_false(intermediate_value_conforms(d3 * (1 + 1e-7), d3))
  # The definitions give the stated k3 within the same two multiples.
  expect_lte(abs(k[["k3"]] / 2^-1074 - k3_intermediate_multiple), 2)
})

# Requirement: Numerical range
# The values that the Example "errors from intermediate values" names. This
# test makes no call of that Example, so it does not carry the Example's
# marker; the marked test is the one above.
test_that("The values named for the calls that stop on an intermediate follow from the definitions", {
  # Where the package stops on both calls of the Example (as it does at the
  # commit that added this test, on R 4.6.1, aarch64 macOS), the marked test
  # above has no returned value to compare. Here each named value is computed
  # apart from the package, from the definitions spec_cumulants() and
  # spec_third_order() above.
  cumulants <- spec_cumulants
  third_order <- spec_third_order
  # The local definitions are checked first: against design_cumulants() on
  # three items whose discriminations are not 1 and whose probabilities at
  # theta are not 1/2 (0.63, 0.55 and 0.17), so that every power of a and every
  # odd cumulant enters; then the formula on two rows of the Example "Delta_3
  # table", against the exact values of the table and against delta3().
  a <- c(0.7, 1.3, 2.1)
  b <- c(-0.4, 0.2, 1.1)
  expect_lt(max(abs(cumulants(a, b, 0.35) / design_cumulants(a, b, 0.35) - 1)), 1e-12)
  rows <- list(
    list(a = c(1, 1), b = c(log(5), -log(5)), theta = 0, exact = 110889 / 625),
    list(a = 1, b = 0, theta = stats::qlogis(3 / 10),
         exact = -369968225 / 583443)
  )
  for (row in rows) {
    local <- third_order(cumulants(row$a, row$b, row$theta))
    expect_lt(abs(local / row$exact - 1), 1e-10)
    expect_lt(abs(local / delta3(row$a, row$b, row$theta, "MLE") - 1), 1e-10)
  }

  k <- cumulants(c(1, 1000), c(730, 0), 0)
  # k3: the item at its difficulty contributes exactly 0; the other item has
  # P = x / (1 + x) with x = exp(-730), so u (1 - 2P) = x (1 - x) / (1 + x)^3.
  # The named value is a subnormal number, where neighbouring doubles are
  # 2^-1074 apart, a relative 1.1e-6 of this value, so it is compared within
  # two such spacings. That the cumulant of such a design is exp(eta) / 2
  # more closely than that spacing shows is measured on the same design with
  # the difficulty 700 in place of 730, where exp(-700) / 2 is a normal number.
  named_k3 <- exp(-730) / 2
  expect_true(named_k3 > 0 && named_k3 < 2^-1022)
  expect_lte(abs(k[["k3"]] - named_k3), 2 * 2^-1074)
  at_700 <- cumulants(c(1, 1000), c(700, 0), 0)[["k3"]]
  expect_lt(abs(at_700 / (exp(-700) / 2) - 1), 1e-12)
  # "about 4.6e-318": within half a unit of the last digit given.
  expect_lt(abs(named_k3 - 4.6e-318), 0.05e-318)
  # Delta_3: with k3, k5 and k7 below 1e-317, only the terms in k4 k6 and
  # k4^3 remain, and the value is -1936 / 3e6.
  expect_lt(abs(third_order(k) / (-1936 / 3e6) - 1), 1e-12)
})

# Whether `out`, a returned scaled risk or the error of the call, is a finite
# risk that is not negative or an error naming `field` as not representable.
finite_risk_or_stop <- function(out, field) {
  if (inherits(out, "error")) {
    message <- conditionMessage(out)
    return(grepl("representable", message, fixed = TRUE) &&
             grepl(field, message, fixed = TRUE))
  }
  length(out) == 1L && is.finite(out) && out >= 0
}

# Requirement: Numerical range
# Scenario: A risk whose reports agree with the target
test_that("A risk whose reports agree with the target may round to zero", {
  value <- risk_diff(c(1, 1 + 2^-36), c(0, 0), 0, "MLE")$smse_boot
  expect_true(is.finite(value) && value >= 0 && value < 1e-100)
  # The ability one double above -r: the exact risk for the supplied numbers is
  # 7.9e3 multiples of 2^-1074, and the reports and the target can still be
  # computed as the same double.
  r <- 0.6761242835860889
  out <- tryCatch(risk_diff(rep(1000, 3), c(-2 * r, 0, 2 * r), -r + 2^-53,
                            "MLE")$smse_boot, error = identity)
  expect_true(finite_risk_or_stop(out, "smse_boot"))
  expect_false(finite_risk_or_stop(-2^-1074, "smse_boot"))
  expect_false(finite_risk_or_stop(Inf, "smse_boot"))
  expect_false(finite_risk_or_stop(simpleError("scaled risk smse_info is not representable"),
                                   "smse_boot"))
  expect_false(finite_risk_or_stop(simpleError("scaled risk smse_boot failed"), "smse_boot"))
  expect_true(finite_risk_or_stop(
    simpleError("scaled risk smse_boot is not representable in double precision."),
    "smse_boot"))
})

# Requirement: Numerical range
# Scenario: A design average whose item polynomial is computed as 0
test_that("An item polynomial computed as zero gives zero or a value below 1e-15", {
  zero_or_tiny <- function(x) length(x) == 1L && is.finite(x) && abs(x) < 1e-15
  expect_true(zero_or_tiny(design_cumulants(1, 0, acosh(2))[["k4"]]))
  expect_true(zero_or_tiny(design_cumulants(1, 0, 2.2924316695611773)[["k5"]]))
  expect_false(zero_or_tiny(2e-15))
  expect_false(zero_or_tiny(NaN))
})

# Requirement: Numerical range
test_that("The internal logit guard rejects unrepresentable logits directly", {
  for (eta in c(Inf, -Inf, NaN, NA_real_, 2^-1074, -2^-1074)) {
    expect_error(irtse:::checked_logits(c(0, eta)), "representable", fixed = TRUE)
  }
  ordinary <- c(-720, -1, 0, 1e-18, 720)
  expect_identical(irtse:::checked_logits(ordinary), ordinary)
})

# Requirement: Numerical range; Diagnostic without choosing
test_that("Diagnostic windows may extend beyond the supported input range", {
  for (location in c(0, 1000)) {
    expect_no_error(result <- se_diagnostic(1e-3, location, location, "WLE"))
    expect_lt(abs(result$se_info / 2000 - 1), 1e-12)
  }
})

# Requirement: Package metadata and documentation; Numerical range
test_that("README and help retain both item-level accuracy qualifications", {
  texts <- list()
  readme <- testthat::test_path("..", "..", "README.md")
  if (file.exists(readme)) texts$README <- readLines(readme, warn = FALSE)
  for (name in c("se_info", "se_boot", "design_cumulants", "delta2", "delta3",
                 "risk_diff", "se_diagnostic", "se_select")) {
    rd_name <- paste0(name, ".Rd")
    rd_path <- testthat::test_path("..", "..", "man", rd_name)
    texts[[name]] <- if (file.exists(rd_path)) readLines(rd_path, warn = FALSE) else {
      as.character(tools::Rd_db("irtse")[[rd_name]])
    }
  }
  for (name in names(texts)) {
    text <- gsub("\\s+", " ", paste(texts[[name]], collapse = " "))
    expect_match(text, "cumulant polynomial near its zeros", fixed = TRUE, info = name)
    expect_match(text, "theta - b", fixed = TRUE, info = name)
  }
})

# The text of a document after markup is removed, as the Example "documents on
# the limits of the double range" describes. White space is the six ASCII
# white-space characters, named one by one so that the reading does not depend
# on the locale. A code span is replaced in one pass, so in a code span inside
# another only the inner one is replaced. A known limit: a single code span of
# about five million characters exceeds the matching limit of the
# regular-expression engine, so that span is not read; each document of the
# package has fewer than 16,000 characters (README.md, the longest, has about
# 15,200).
normalise_document_text <- function(text) {
  text <- paste(text, collapse = " ")
  text <- gsub("(?s)\\\\code\\{((?:[^{}\\\\]|\\\\.)*?)\\}", "\\1", text, perl = TRUE)
  text <- gsub("\\\\([%{}^])", "\\1", text, perl = TRUE)
  text <- gsub("`", "", text, fixed = TRUE)
  gsub("[ \t\n\v\f\r]+", " ", text)
}

# A sentence of a document, checked as fixed text: the text from its first word
# to its full stop has to stand at the start of the document text or directly
# after a full stop and a space. This is a check of text, not of meaning: the
# same text matches after any full stop and a space (an abbreviation, a
# numbered item, a sentence that negates it, a comment line of a help page read
# from its source file) and with more text directly after its full stop.
has_document_sentence <- function(text, sentence) {
  at <- gregexpr(sentence, text, fixed = TRUE)[[1]]
  if (at[1] < 0) return(FALSE)
  any(vapply(at, function(i) i == 1 || identical(substr(text, i - 2, i - 1), ". "),
             logical(1)))
}

# The sentences that README.md and the help pages of the eight exported
# functions state, as the Example "documents on the limits of the double
# range" lists them: the text after markup is removed.
risk_limit_sentences <- c(
  terms = paste0(
    "The four scaled risks (smse_info, smse_boot, diff, smse_select) are summed",
    " on the log scale: each term is formed as a whole from the logarithm of ",
    "its weight and of its two differences, so the sum never forms a weight, a ",
    "squared difference or a partial product on its own, and none of them can ",
    "stop a call or remove a term by leaving the double range."),
  checked = paste0(
    "Apart from requiring that the log weights and differences it receives are ",
    "finite, the sum checks the range only on the logarithm of the whole sum ",
    "and on the value rebuilt from it."),
  stop = paste0(
    "A risk whose logarithm, as this summation computes it, is below ",
    "log(2^-1074) stops the call with representable and the field name, even ",
    "when ordinary arithmetic would round the value up to 2^-1074; a risk whose",
    " value rebuilt from that logarithm is not finite, or is 0, stops the call ",
    "in the same way."),
  rebuilt = paste0(
    "The value of a scaled risk is rebuilt from its logarithm."),
  attribution = paste0(
    "Rounding that logarithm alone contributes a relative error of the order of",
    " 1e-13 where the logarithm is near 700 in absolute value; errors already ",
    "in the reports and log weights add to it, and so does rounding the result ",
    "to a double."),
  subnormal = paste0(
    "A returned risk in the subnormal range is a multiple of 2^-1074, so ",
    "rounding the result to a subnormal double can give a much larger relative ",
    "error: in risk_diff(rep(1000, 3), c(-1.44, 0, 1.44), -0.72 + 5e-9, \"MLE\") ",
    "the exact smse_boot is 1.6 times 2^-1074 and the call returned 2 times ",
    "2^-1074."),
  no_bound = paste0(
    "A subnormal risk is not claimed to be within 2^-1074 of the exact value, ",
    "and no general bound on the relative error of a returned risk is claimed."),
  follow = paste0(
    "Both decisions follow the logarithm computed from the reports, the target ",
    "and the log weights in double precision, not the exact risk, and the two ",
    "can lie on different sides of 2^-1074 or of the largest double."),
  not_close = paste0(
    "At 2^-1074 this is not confined to risks close to the limit: when a report",
    " and the target nearly agree, their difference keeps few significant ",
    "digits and the computed risk can be wrong by a large factor."),
  call_stopped = paste0(
    "In risk_diff(rep(1000, 3), c(-2 * r, 0, 2 * r), -r + 2^-52, \"MLE\") with r ",
    "set to 0.685794632511088 the report and the target agree to a relative ",
    "1.1e-13, the exact smse_boot for the supplied numbers is 2.05 times ",
    "2^-1074, which is representable, and the call stopped."),
  call_returned = paste0(
    "With r set to 0.6907256300571373 and the ability -r + 2^-49 the exact ",
    "value for the supplied numbers is 0.96 times 2^-1074, which is below the ",
    "smallest positive double, and the call returned 2^-1074."),
  upper = paste0(
    "At the largest double, in risk_diff(rep(512, 2), c(-B, B), 0, \"MLE\") with ",
    "B set to 1.4106628166864512, the double that R gives for 1042 * log(2) / ",
    "512, the exact smse_info is above the largest double by a relative ",
    "9.9e-15 and the call returned 1.7976931348622732e308."),
  no_claim = paste0(
    "In general the package does not claim that the outcome agrees with the ",
    "side of either limit on which the exact risk lies."),
  platform = paste0(
    "The outcomes stated for the calls at the abilities -0.72 + 5e-9, -r + ",
    "2^-52 and -r + 2^-49 and for the call with B were measured with R 4.6.1 ",
    "on aarch64 macOS and can differ on another platform."),
  zero_diff = paste0(
    "A diff of 0 means the two risks are not distinguished in double precision,",
    " not that they are equal."),
  diff_stop = paste0(
    "A nonzero diff whose logarithm of the absolute value, computed in the same",
    " way, is below log(2^-1074) stops the call with representable and the ",
    "field name."))

# The attribution of the error of the order of 1e-13, through the words after
# the number.
risk_log_rounding_clause <- paste0(
  "Rounding that logarithm alone contributes a relative error of the order of ",
  "1e-13 where the logarithm is near 700 in absolute value;")

# Strings of earlier versions that the documents no longer contain.
retired_risk_sentences <- c(
  "near the ends of the double range its relative error is of the order of 1e-13",
  "within a relative 1e-12 of 2^-1074",
  "inside the band at the two ends of the double range",
  "A positive risk below 2^-1074 stops the call",
  "A nonzero diff with absolute value below 2^-1074 stops the call",
  "A risk that is itself a subnormal number",
  "a risk of 1.6 such spacings is returned as 2",
  "This is not confined to risks close to those limits",
  "Beyond the calls named in the specification",
  "a value above the largest double stops the call",
  "1.6 such spacings",
  "the exact smse_boot is 2.05 times 2^-1074",
  "the exact value is 0.96 times 2^-1074",
  "with B set to 1042 * log(2) / 512 the exact",
  "The outcomes of the calls named here")

# Whether the tests run in the source folder of the package: the folder two
# levels up is a source tree (is_source_tree(): an R folder and no file
# Meta/package.rds), as opposed to an installed package, also one installed with its
# tests, or the folder of R CMD check. A document that the source folder
# should have is then required, not skipped.
in_source_folder <- function() is_source_tree(testthat::test_path("..", ".."))

# Requirement: Numerical range
# A control of the predicate that the Example "documents on the limits of the
# double range" states for the source folder; the specification has no
# Scenario of this name.
test_that("The source folder is told apart from an installed package", {
  base <- tempfile("irtse-tree-")
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  tree <- function(name, folders, files = character()) {
    path <- file.path(base, name)
    dir.create(path, recursive = TRUE)
    for (folder in folders) dir.create(file.path(path, folder))
    for (file in files) {
      dir.create(dirname(file.path(path, file)), recursive = TRUE, showWarnings = FALSE)
      writeLines("A stand-in.", file.path(path, file))
    }
    path
  }
  expect_true(is_source_tree(tree("source", "R")))
  # A Meta folder without Meta/package.rds: devtools writes a vignette index
  # into one in a source folder with vignettes, and on a case-insensitive file
  # system a folder named meta also answers to Meta.
  expect_true(is_source_tree(tree("vignette_index", "R", "Meta/vignette.rds")))
  expect_true(is_source_tree(tree("lower_case_meta", c("R", "meta"))))
  expect_false(is_source_tree(tree("installed", "R", "Meta/package.rds")))
  expect_false(is_source_tree(tree("package_rds_only", character(), "Meta/package.rds")))
  expect_false(is_source_tree(tree("neither", character())))
})

# The exported functions whose help pages carry the template of the numerical
# limits.
template_pages <- c("se_info", "se_boot", "design_cumulants", "delta2", "delta3",
                    "risk_diff", "se_diagnostic", "se_select")

# A help page as text: its source file where the source folder is present, or,
# in an installed package, the page from the Rd database `db`. The pieces of
# the database are joined without a separator; a space between them would
# stand between a code span and the comma, semicolon or bracket next to it.
help_page_text <- function(name, db) {
  if (in_source_folder()) {
    rd_path <- testthat::test_path("..", "..", "man", paste0(name, ".Rd"))
    return(readLines(rd_path, warn = FALSE))
  }
  paste(as.character(db[[paste0(name, ".Rd")]]), collapse = "")
}

# The rows of the Example "documents on the limits of the double range" for one
# document, given as its lines.
expect_risk_limit_document <- function(lines, name) {
  text <- gsub("[ \t\n\v\f\r]+", " ", paste(lines, collapse = " "))
  expect_match(text, "diff[^.]*not distinguished[^.]*double precision", info = name)
  text <- normalise_document_text(lines)
  for (sentence in names(risk_limit_sentences)) {
    expect_true(has_document_sentence(text, risk_limit_sentences[[sentence]]),
                info = paste(name, sentence))
  }
  expect_match(text, risk_log_rounding_clause, fixed = TRUE, info = name)
  expect_match(text, paste0("rounding the result to a subnormal double can give ",
                            "a much larger relative error"), fixed = TRUE, info = name)
  expect_match(text, "A returned risk in the subnormal range is a multiple of 2^-1074",
               fixed = TRUE, info = name)
  for (retired in retired_risk_sentences) {
    expect_false(grepl(retired, text, fixed = TRUE), info = paste(name, retired))
  }
}

# Requirement: Numerical range
# Scenario: Documents state the rules at the limits of the double range
# spec-example: irtse-standard-errors / documents on the limits of the double range
test_that("README states the rules at the limits of the double range", {
  if (!in_source_folder()) skip("README.md is checked only in the source folder of the package")
  readme <- testthat::test_path("..", "..", "README.md")
  expect_true(file.exists(readme))
  expect_risk_limit_document(readLines(readme, warn = FALSE), "README")
})

# Requirement: Numerical range
# Scenario: Documents state the rules at the limits of the double range
# spec-example: irtse-standard-errors / documents on the limits of the double range
test_that("The eight help pages state the rules at the limits of the double range", {
  db <- if (in_source_folder()) NULL else tools::Rd_db("irtse")
  for (name in template_pages) expect_risk_limit_document(help_page_text(name, db), name)
})

# Requirement: Numerical range
# Scenario: Documents state the rules at the limits of the double range
test_that("The document sentence check compares fixed text from first word to full stop and reads through markup", {
  wanted <- risk_limit_sentences[["follow"]]
  rd <- paste0(
    "Both decisions follow the logarithm computed from the reports, the target and ",
    "the log weights in double precision, not the exact risk, and the two can lie ",
    "on different sides of \\code{2^-1074} or of the largest double.")
  markdown <- sub("\\code{2^-1074}", "`2^-1074`", rd, fixed = TRUE)
  check <- function(raw) has_document_sentence(normalise_document_text(raw), wanted)
  expect_true(check(rd))
  expect_true(check(markdown))
  expect_true(check(paste("A sentence before it.", rd, "A sentence after it.")))
  expect_false(check(paste("It is false that", rd)))
  expect_false(check(sub("largest double.", "largest double only for diff.", rd,
                         fixed = TRUE)))
  expect_false(check(sub(", not the exact risk", "", rd, fixed = TRUE)))
  expect_false(check(sub("2^-1074", "12^-1074", rd, fixed = TRUE)))

  # What a check of fixed text does not see, as the Example states it: the
  # sentence after any full stop and a space (an abbreviation, a numbered item,
  # a sentence that negates it, text in which the full stop and the space meet
  # only once backticks are dropped, a comment line of a help page read from
  # its source file), and with more text directly after its full stop. An HTML
  # comment, which the Example "document writing restrictions" keeps out of
  # README.md, is read the same way.
  expect_true(check(paste("This is false, e.g.", rd)))
  expect_true(check(paste("1.", rd)))
  expect_true(check(paste("None of this holds.", rd)))
  expect_true(check(paste0("Retracted.` `", rd)))
  expect_true(check(c("% An earlier text.", rd)))
  expect_true(check(paste0("<!-- An earlier text. ", markdown, " -->")))
  expect_true(check(sub("largest double.", "largest double.5 only", rd, fixed = TRUE)))

  # White space is the six ASCII white-space characters, whatever the locale:
  # any of them between two words reads as one space, and no other character
  # does, so a sentence written with a no-break, thin or ideographic space is
  # not found.
  spaced <- function(ch) sub("the log weights", paste0("the", ch, "log weights"), rd,
                             fixed = TRUE)
  for (ch in c("\t", "\n", "\v", "\f", "\r", "  ")) {
    expect_true(check(spaced(ch)), info = sprintf("U+%04X", utf8ToInt(substr(ch, 1, 1))))
  }
  for (ch in c("\u00a0", "\u2009", "\u3000")) {
    expect_false(check(spaced(ch)), info = sprintf("U+%04X", utf8ToInt(ch)))
  }

  # A code span is replaced in one pass: in a code span inside another only the
  # inner one is replaced and the outer markup stays in the text.
  expect_identical(normalise_document_text("A \\code{\\code{2^-1074}} B"),
                   "A \\code{2^-1074} B")

  # A code span is replaced by its content as it stands: white space inside the
  # braces is kept on both sides, so a space that splits a word, or that
  # completes a retired sentence, is seen.
  expect_false(check(sub("largest double", "lar\\code{ gest} double", rd, fixed = TRUE)))
  expect_true(grepl(
    "A positive risk below 2^-1074 stops the call",
    normalise_document_text("A positive risk below\\code{ 2^-1074} stops the call."),
    fixed = TRUE))
  expect_true(grepl(
    "A positive risk below 2^-1074 stops the call",
    normalise_document_text("A positive risk below \\code{2^-1074 }stops the call."),
    fixed = TRUE))
  # The content may hold a backslash directly in front of a line break. The
  # text of an installed help page keeps its line breaks, so the span has to be
  # replaced there as it is where the lines of a source file are joined.
  expect_true(grepl(
    "A positive risk below 2^-1074 stops the call",
    normalise_document_text("A positive risk below\\code{ 2^-1074 stops the call\\\n} is here."),
    fixed = TRUE))
  # The same span given as two lines of a source file, as readLines() returns
  # them.
  expect_true(grepl(
    "A positive risk below 2^-1074 stops the call",
    normalise_document_text(c("A positive risk below\\code{ 2^-1074 stops the call\\",
                              "} is here.")),
    fixed = TRUE))
  # With white space between the macro and its brace R reads no code span and
  # reports the page; when a space or a tab stands directly in front of the
  # brace it also drops the rest of the section. The text stays as written and
  # completes no sentence.
  expect_false(check(sub("\\code{2^-1074}", "\\code {2^-1074}", rd, fixed = TRUE)))

  clause <- function(raw) {
    grepl(risk_log_rounding_clause, normalise_document_text(raw), fixed = TRUE)
  }
  written <- paste0(
    "Rounding that logarithm alone contributes a relative error of the order of ",
    "1e-13 where the logarithm is near 700 in absolute value; errors already in")
  expect_true(clause(written))
  expect_false(clause(sub("1e-13 ", "1e-130 ", written, fixed = TRUE)))
  expect_false(clause(sub("1e-13 ", "1e-13.5 ", written, fixed = TRUE)))
  expect_false(clause(sub("near 700 in absolute value;", "near 700;", written,
                          fixed = TRUE)))
})

# Requirement: Numerical range
# Scenario: Documents state the rules at the limits of the double range
test_that("No source help page contains a retired sentence on risks near the limits", {
  man <- testthat::test_path("..", "..", "man")
  if (!dir.exists(man)) skip("Source help pages not available in installed package")
  paths <- list.files(man, pattern = "[.]Rd$", full.names = TRUE)
  bad <- vapply(paths, function(path) {
    text <- normalise_document_text(readLines(path, warn = FALSE))
    any(vapply(retired_risk_sentences, function(retired) {
      grepl(retired, text, fixed = TRUE)
    }, logical(1)))
  }, logical(1))
  expect_false(any(bad), info = paste(basename(paths[bad]), collapse = ", "))
})

# The restrictions of the Example "document writing restrictions" that one
# document breaks. `kind` is "help" for a help page or the help template, and
# "readme" for README.md. The ASCII row tests `bytes`, the bytes of the file;
# the other rows test `lines`, the lines that readLines() returns, which are
# what the sentence check reads (readLines() ends a line at a NUL). For a
# planted text given only as lines, the bytes are those of the lines.
document_restriction_failures <- function(lines, kind,
                                          bytes = charToRaw(paste(lines, collapse = "\n"))) {
  failures <- character()
  if (any(as.integer(bytes) > 127L)) {
    failures <- c(failures, "not ASCII")
  }
  if (kind == "help") {
    if (any(grepl("\\\\[\\\\{}]", lines, perl = TRUE))) {
      failures <- c(failures, "backslash before a backslash or a brace")
    }
    if (any(grepl("\\\\$", lines, perl = TRUE))) {
      failures <- c(failures, "backslash at the end of a line")
    }
    if (any(grepl("\\\\code(?!\\{)", lines, perl = TRUE))) {
      failures <- c(failures, "code macro not followed by a brace")
    }
  } else {
    if (any(grepl("<!--", lines, fixed = TRUE))) failures <- c(failures, "HTML comment")
    if (any(grepl("~~", lines, fixed = TRUE))) failures <- c(failures, "strikethrough")
  }
  failures
}

# Requirement: Numerical range
# Scenario: Documents are written in the form that the sentence check reads
# spec-example: irtse-standard-errors / document writing restrictions
test_that("The documents keep out the forms that the sentence check does not read as R does", {
  # Each restriction reports a planted text that breaks it, and only that one.
  expect_identical(document_restriction_failures("caf\u00e9", "readme"), "not ASCII")
  expect_identical(document_restriction_failures("caf\u00e9", "help"), "not ASCII")
  for (text in c("a \\\\ b", "a \\{ b", "a \\} b")) {
    expect_identical(document_restriction_failures(text, "help"),
                     "backslash before a backslash or a brace", info = text)
  }
  expect_identical(document_restriction_failures(c("a b\\", "c"), "help"),
                   "backslash at the end of a line")
  for (text in c("\\code {x}", "\\code", "\\code\tx")) {
    expect_identical(document_restriction_failures(text, "help"),
                     "code macro not followed by a brace", info = text)
  }
  expect_identical(document_restriction_failures("<!-- a -->", "readme"), "HTML comment")
  expect_identical(document_restriction_failures("~~a~~", "readme"), "strikethrough")
  expect_identical(document_restriction_failures("\\code{x} and `y`, 50\\%.", "help"),
                   character())
  # The ASCII row reads the bytes of a file: a byte-order mark, and bytes after
  # a NUL on the same line, where readLines() ends the line, are seen.
  bom <- tempfile(fileext = ".md")
  nul <- tempfile(fileext = ".md")
  on.exit(unlink(c(bom, nul)), add = TRUE)
  writeBin(c(as.raw(c(0xef, 0xbb, 0xbf)), charToRaw("Plain text.\n")), bom)
  writeBin(c(charToRaw("A"), as.raw(0), charToRaw(" caf"), as.raw(c(0xc3, 0xa9)),
             charToRaw("\n")), nul)
  restriction_failures_of_file <- function(path, kind) {
    document_restriction_failures(readLines(path, warn = FALSE), kind,
                                  readBin(path, "raw", file.info(path)$size))
  }
  for (path in c(bom, nul)) {
    expect_identical(restriction_failures_of_file(path, "readme"), "not ASCII",
                     info = basename(path))
  }
  # The other rows read the lines that readLines() returns, and a line ends at
  # a NUL: an HTML comment that follows a NUL on its line is not seen by them.
  after_nul <- tempfile(fileext = ".md")
  on.exit(unlink(after_nul), add = TRUE)
  writeBin(c(charToRaw("A"), as.raw(0), charToRaw(" <!-- a -->\n")), after_nul)
  expect_identical(restriction_failures_of_file(after_nul, "readme"), character())

  if (!in_source_folder()) {
    skip("The source files of the documents are checked only in the source folder of the package")
  }
  root <- testthat::test_path("..", "..")
  template <- file.path(root, "man-roxygen", "numerical-accuracy.R")
  readme <- file.path(root, "README.md")
  expect_true(file.exists(template))
  expect_true(file.exists(readme))
  pages <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
  expect_length(pages, 9L)
  for (path in c(pages, template)) {
    expect_identical(restriction_failures_of_file(path, "help"), character(),
                     info = basename(path))
  }
  expect_identical(restriction_failures_of_file(readme, "readme"), character())
})

# TRUE on the version and platform on which the outcomes of the Example
# "documented outcomes" were measured.
documented_platform <- function(major = R.version$major, minor = R.version$minor,
                                arch = R.version$arch,
                                sysname = Sys.info()[["sysname"]]) {
  identical(c(major, minor, arch, sysname), c("4", "6.1", "aarch64", "Darwin"))
}

# Requirement: Numerical range
# Scenario: The documented outcomes on the platform they were measured on
# spec-example: irtse-standard-errors / documented outcomes
test_that("The documented outcomes hold with R 4.6.1 on aarch64 macOS", {
  expect_true(documented_platform("4", "6.1", "aarch64", "Darwin"))
  expect_false(documented_platform("5", "6.1", "aarch64", "Darwin"))
  expect_false(documented_platform("4", "6.2", "aarch64", "Darwin"))
  expect_false(documented_platform("4", "6.1", "x86_64", "Darwin"))
  expect_false(documented_platform("4", "6.1", "aarch64", "Linux"))
  if (!documented_platform()) {
    skip("The documented outcomes were measured with R 4.6.1 on aarch64 macOS")
  }
  q <- 2^-1074
  expect_identical(risk_diff(rep(1000, 3), c(-1.44, 0, 1.44), -0.72 + 5e-9,
                             "MLE")$smse_boot, 2 * q)
  r <- 0.685794632511088
  stopped <- tryCatch(risk_diff(rep(1000, 3), c(-2 * r, 0, 2 * r), -r + 2^-52, "MLE"),
                      error = identity)
  expect_s3_class(stopped, "error")
  expect_match(conditionMessage(stopped), "representable", fixed = TRUE)
  expect_match(conditionMessage(stopped), "smse_boot", fixed = TRUE)
  r <- 0.6907256300571373
  expect_identical(risk_diff(rep(1000, 3), c(-2 * r, 0, 2 * r), -r + 2^-49,
                             "MLE")$smse_boot, q)
  # 1.7976931348622732e308 to 17 significant digits, written in hexadecimal:
  # R does not round every decimal literal correctly, and reads that one as a
  # double six spacings higher.
  B <- 1042 * log(2) / 512
  expect_identical(risk_diff(rep(512, 2), c(-B, B), 0, "MLE")$smse_info,
                   0x1.fffffffffff2ap+1023)
  r <- 0.6761242835860889
  expect_identical(risk_diff(rep(1000, 3), c(-2 * r, 0, 2 * r), -r + 2^-53,
                             "MLE")$smse_boot, 0)
  expect_identical(design_cumulants(1, 0, acosh(2))[["k4"]], 0)
  expect_identical(design_cumulants(1, 0, 2.2924316695611773)[["k5"]], 0)
})

# Requirement: Numerical range
# Scenario: The outcome at the limits of the double range follows the computed logarithm
# spec-example: irtse-standard-errors / outcomes at the smallest positive double
test_that("Calls with an exact risk near the smallest positive double stop or return a finite risk", {
  calls <- list(
    list(r = 0.685794632511088, theta = -0.685794632511088 + 2^-52),
    list(r = 0.6907256300571373, theta = -0.6907256300571373 + 2^-49),
    list(r = 0.7204421116605104, theta = -0.7204421067316231),
    list(r = 0.7198038989710372, theta = -0.7198038953855531))
  for (call in calls) {
    out <- tryCatch(
      risk_diff(rep(1000, 3), c(-2 * call$r, 0, 2 * call$r), call$theta, "MLE"),
      error = function(e) e)
    label <- format(call$r, digits = 17)
    if (inherits(out, "error")) {
      expect_match(conditionMessage(out), "representable", fixed = TRUE, info = label)
      expect_match(conditionMessage(out), "smse_boot", fixed = TRUE, info = label)
    } else {
      expect_true(is.finite(out$smse_boot), info = label)
      expect_gte(out$smse_boot, 0)
    }
  }
})

# Requirement: Package metadata and documentation
test_that("The scan for local paths reads package text files only", {
  paths <- c("pkg/R/irtse.rdb", "pkg/R/irtse.rdx", "pkg/Meta/package.rds",
             "pkg/help/aliases.RDS", "pkg/.git", "pkg/.git/logs/HEAD",
             "pkg/.git/hooks/pre-commit.sample", ".git/config",
             "pkg/R/irtse", "pkg/NEWS.md", "pkg/R/risk.R", "pkg/man/se_boot.Rd",
             "pkg/.Rbuildignore", "pkg/tests/testthat/words.rds.txt",
             "pkg/R/.gitkeep")
  expect_identical(scannable_files(paths),
                   c("pkg/R/irtse", "pkg/NEWS.md", "pkg/R/risk.R",
                     "pkg/man/se_boot.Rd", "pkg/.Rbuildignore",
                     "pkg/tests/testthat/words.rds.txt", "pkg/R/.gitkeep"))
  # Three bytes of a compressed index read as a drive path: the chance match
  # that the scan must not meet.
  compressed <- as.raw(c(0x34, 0x6e, 0x45, 0x3a, 0x5c, 0x58))
  drive_path <- paste0("[A-Za-z]", ":", "\\\\")
  expect_gt(length(grepRaw(drive_path, compressed)), 0L)
})

# Requirement: Package metadata and documentation
test_that("POSIX local paths start at a path boundary and exclude URL paths", {
  roots <- c("Users", "home", "tmp", "var", "private", "opt", "Volumes",
             "srv", "etc", "mnt", "usr", "root")
  matches <- function(text) {
    any(vapply(additional_local_path_patterns(), function(pattern) {
      length(grepRaw(pattern, charToRaw(text))) > 0L
    }, logical(1)))
  }
  for (root in roots) {
    path <- paste0("/", root, "/example.txt")
    for (prefix in c("", " ", "\n", "'", "(", "=")) {
      expect_true(matches(paste0(prefix, path)), info = paste(root, prefix))
    }
    for (prefix in c("a", "Z", "1", ".", "-", "_", "/", ":",
                     "https://example.org")) {
      expect_false(matches(paste0(prefix, path)), info = paste(root, prefix))
    }
  }
})

# Requirement: Numerical range
test_that("A positive risk below the smallest positive double is rejected", {
  # 3 r^2 (sech(s r / 2) - sech(s theta / 2))^2 has log risk
  # below log(2^-1074), although the reports differ from the target.
  expect_error(risk_diff(rep(1000, 3), c(-1.44, 0, 1.44),
                         -0.72 + 1e-9, "MLE"), "representable", fixed = TRUE)
  expect_error(risk_diff(rep(512, 3), c(-2.8125, 0, 2.8125),
                         -1.40625 + 2^-30, "MLE"), "representable", fixed = TRUE)
  rounded <- risk_diff(c(1, 1 + 2^-36), c(0, 0), 0, "MLE")$smse_boot
  expect_true(is.finite(rounded) && rounded >= 0 && rounded < 1e-100)
})

# Requirement: Package metadata and documentation
test_that("Source scans retain serialized suffixes but installed scans omit them", {
  paths <- c("pkg/tests/words.rds", "pkg/R/words.rdb", "pkg/R/words.rdx")
  source <- tryCatch(scannable_files(paths, installed = FALSE),
                     error = conditionMessage)
  expect_identical(source, paths)
  expect_identical(scannable_files(paths), character())
})

# Requirement: Package metadata and documentation
test_that("Local path scans include unknown roots and exclude URL paths", {
  path <- paste0("/", "unlisted-root", "/private.txt")
  matches <- function(text) {
    any(vapply(additional_local_path_patterns(), function(pattern) {
      length(grepRaw(pattern, charToRaw(text))) > 0L
    }, logical(1)))
  }
  expect_true(matches(path))
  expect_true(matches(paste0("path=", path)))
  expect_false(matches(paste0("https://example.org", path)))
  expect_false(matches(paste0("relative", path)))
})

# Requirement: Package metadata and documentation (package files contain no absolute local paths)
# Design decision of the change irtse-risk-accuracy-limits: the path pattern accepts spaces inside the first
# folder name. This is a test of the scan itself; the specification has no Scenario of this name.
test_that("Local path scans detect first folders with spaces without matching prose", {
  matches <- function(text) {
    any(vapply(additional_local_path_patterns(), function(pattern) {
      length(grepRaw(pattern, charToRaw(text))) > 0L
    }, logical(1)))
  }
  for (folder in c(paste("Project", "Files"), paste("Shared", "Project", "Files"),
                   paste("Shared", "Research", "Project", "Data", "Files"),
                   paste("Project", "", "Files"), paste("2", "Project", "Files"),
                   paste("Project", "&", "Files"))) {
    expect_true(matches(paste0("/", folder, "/", "private.txt")))
  }
  controls <- c(paste("a", "/", "b", "/", "c"),
                paste0("1", "/", "2 and 3", "/", "4"),
                paste("k_r", "/", "k_2"),
                paste0("+", "/", "- 1 and x", "/", "y"),
                paste0("https:", "/", "/", "example.org", "/",
                       paste("Project", "Files"), "/", "private.txt"),
                # Two closing HTML tags with text between them, as in an installed help index.
                paste0("<", "/", "div><h2>Documentation for package irtse<", "/", "h2>"))
  for (text in controls) expect_false(matches(text))
})

# Requirement: Package metadata and documentation
test_that("R CMD build excludes a submodule git file", {
  if (!in_source_folder()) {
    skip("The build ignore file is checked only in the source folder of the package")
  }
  path <- testthat::test_path("..", "..", ".Rbuildignore")
  expect_true(file.exists(path))
  patterns <- readLines(path, warn = FALSE)
  expect_true(any(vapply(patterns, function(pattern) {
    grepl(pattern, ".git", perl = TRUE, ignore.case = TRUE)
  }, logical(1))))
})

# Requirement: Package metadata and documentation
# A control of the two scans themselves; the specification has no Scenario of
# this name.
test_that("The scans skip the files that testthat leaves after a failing run", {
  base <- tempfile("irtse-leftover-")
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  put <- function(relative, text) {
    path <- file.path(base, relative)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(text, path)
  }
  local_path <- paste0("/", "Users", "/", "someone", "/", "notes.txt")
  put("DESCRIPTION", "Package: control")
  put("README.md", "A control tree.")
  put("NEWS.md", "No news.")
  put(file.path("R", "code.R"), "f <- function() 1")
  put(file.path("tests", "testthat", "testthat-problems.rds"), local_path)
  put(file.path("tests", "testthat", "_problems", "test-f.R"), local_path)
  put(file.path("tests", "testthat", "test-f.R"), paste("# see", local_path))
  # Files of the same names below a folder that ships are still read.
  put(file.path("inst", "tests", "testthat", "_problems", "test-g.R"), local_path)
  put(file.path("inst", "tests", "testthat", "testthat-problems.rds"), local_path)
  planted <- file.path(base, "tests", "testthat", "test-f.R")
  shipped <- file.path(base, "inst", "tests", "testthat",
                       c("_problems/test-g.R", "testthat-problems.rds"))
  expect_setequal(files_matching(shipped_scan_files(base, base), shipped_scan_patterns()),
                  c(planted, shipped))
  # The source scan lists DESCRIPTION, README.md, NEWS.md and the folders R,
  # man and tests; inst is not among them.
  expect_identical(files_matching(source_scan_files(base, base), source_scan_patterns()),
                   planted)
})

# Requirement: Package metadata and documentation
# A control of the two scans themselves; the specification has no Scenario of
# this name.
test_that("The scans read the testthat leftover folder in an installed package", {
  base <- tempfile("irtse-installed-")
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  put <- function(relative, text) {
    path <- file.path(base, relative)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(text, path)
  }
  local_path <- paste0("/", "Users", "/", "someone", "/", "notes.txt")
  # A package installed with its tests has an R folder, as a source folder
  # has, and the file Meta/package.rds, which a source folder does not have.
  # Its folder tests/testthat holds the tests folder of the source (R CMD
  # INSTALL --install-tests), leftovers of a failing run included, or
  # inst/tests/testthat, so a file there ships, whatever its name.
  put("DESCRIPTION", "Package: control")
  put(file.path("Meta", "package.rds"), "A stand-in.")
  put(file.path("R", "control"), "A stand-in.")
  put(file.path("tests", "testthat", "test-f.R"), paste("# see", local_path))
  put(file.path("tests", "testthat", "_problems", "test-g.R"), local_path)
  planted <- file.path(base, "tests", "testthat", c("test-f.R", "_problems/test-g.R"))
  expect_setequal(files_matching(shipped_scan_files(base, base), shipped_scan_patterns()),
                  planted)
  expect_setequal(files_matching(source_scan_files(base, base), source_scan_patterns()),
                  planted)
})

# Requirement: Package metadata and documentation
test_that("R CMD build excludes the files that testthat leaves after a failing run", {
  if (!in_source_folder()) {
    skip("The build ignore file is checked only in the source folder of the package")
  }
  path <- testthat::test_path("..", "..", ".Rbuildignore")
  expect_true(file.exists(path))
  patterns <- readLines(path, warn = FALSE)
  ignored <- function(file) {
    any(vapply(patterns, function(pattern) {
      grepl(pattern, file, perl = TRUE, ignore.case = TRUE)
    }, logical(1)))
  }
  expect_true(ignored("tests/testthat/_problems"))
  expect_true(ignored("tests/testthat/testthat-problems.rds"))
  expect_false(ignored("tests/testthat/test-standard-errors.R"))
  expect_false(ignored("tests/testthat/helper-local-paths.R"))
})

# Requirement: Package metadata and documentation
# The forms that the comment of the path patterns lists as not found; the
# specification has no Scenario of this name.
test_that("Path forms that the scan does not find", {
  matches <- function(text) {
    any(vapply(additional_local_path_patterns(), function(pattern) {
      length(grepRaw(pattern, charToRaw(text))) > 0L
    }, logical(1)))
  }
  s <- "/"
  not_found <- c(
    file_url = paste0("file:", s, s, s, "Users", s, "someone", s, "notes.txt"),
    include_flag = paste0("-I", s, "usr", s, "include"),
    library_flag = paste0("-L", s, "opt", s, "lib"),
    path_entry = paste0("bin:", s, "opt", s, "tools", s, "bin"),
    home_one = paste0("~", s, "notes.txt"),
    home_user = paste0("~someone", s, "notes", s, "a.txt"),
    unlisted_space = paste0(s, "My Folder", s, " sub", s, "file.R"))
  for (name in names(not_found)) expect_false(matches(not_found[[name]]), info = name)
  # A first folder from the list of roots is found even when the second
  # component begins with a space.
  expect_true(matches(paste0(s, "opt", s, " sub", s, "file.R")))
})
