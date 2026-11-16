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
  authors <- eval(parse(text = description$`Authors@R`))
  expect_length(authors, 1L)
  expect_identical(authors$given, "Che")
  expect_identical(authors$family, "Cheng")
  expect_identical(authors$email, "che830621@stat.sinica.edu.tw")
  expect_identical(authors$role, c("aut", "cre"))
})

# Requirement: Package metadata and documentation
test_that("Package files contain no local paths or private issue references", {
  paths <- list.files(system.file(package = "irtse"), recursive = TRUE,
                      full.names = TRUE, all.files = TRUE)
  source_r <- testthat::test_path("..", "..", "R")
  if (dir.exists(source_r)) paths <- c(paths, list.files(source_r, full.names = TRUE))
  paths <- scannable_files(paths, installed = !dir.exists(source_r))
  forbidden_path <- posix_local_path_pattern()
  private_issue <- paste0("Psych", "Quant/[^\r\n]*", "#", "[0-9]+")
  bad <- vapply(unique(paths), function(path) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    length(grepRaw(forbidden_path, bytes)) > 0L ||
      length(grepRaw(private_issue, bytes)) > 0L
  }, logical(1))
  expect_false(any(bad), info = paste(names(bad)[bad], collapse = ", "))
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
  source <- testthat::test_path("..", "..")
  if (dir.exists(file.path(source, "R"))) {
    paths <- c(file.path(source, c("DESCRIPTION", "README.md", "NEWS.md")),
               unlist(lapply(c("R", "man", "tests"), function(dir) {
                 list.files(file.path(source, dir), recursive = TRUE,
                            full.names = TRUE, all.files = TRUE)
               }), use.names = FALSE))
  } else {
    paths <- list.files(system.file(package = "irtse"), recursive = TRUE,
                        full.names = TRUE, all.files = TRUE)
  }
  paths <- scannable_files(paths, installed = !dir.exists(file.path(source, "R")))
  patterns <- c(paste0("[A-Za-z]", ":", "\\\\"),
                paste0("Psych", "Quant/[^\r\n]*", "#", "[0-9]+"),
                paste0("/", "issues", "/", "[0-9]+"))
  patterns <- c(patterns, additional_local_path_patterns())
  bad <- vapply(paths, function(path) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    any(vapply(patterns, function(pattern) length(grepRaw(pattern, bytes)) > 0L,
               logical(1)))
  }, logical(1))
  expect_false(any(bad), info = paste(paths[bad], collapse = ", "))
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

# Requirement: Numerical range
# spec-example: irtse-standard-errors / errors from intermediate values
test_that("Errors from intermediate values may stop an otherwise representable result", {
  cases <- list(
    list(call = function() {
      risk_diff(rep(0.7092, 3), c(-1000, 1000, 1000), 1000, "MLE")$smse_info
    }, expected = exp(709.2) / (2 * sqrt(2) * 0.7092^2)),
    list(call = function() design_cumulants(c(1, 1000), c(730, 0), 0)[["k3"]],
         expected = exp(-730) / 2),
    list(call = function() delta3(c(1, 1000), c(730, 0), 0, "MLE"),
         expected = -1936 / 3e6)
  )
  for (case in cases) {
    value <- tryCatch(case$call(), error = identity)
    if (inherits(value, "error")) {
      expect_match(conditionMessage(value), "representable", fixed = TRUE)
    } else {
      expect_true(length(value) == 1L && is.finite(value) && value != 0)
      expect_lt(abs(value / case$expected - 1), 1e-8)
    }
  }
})

# Requirement: Numerical range
test_that("A risk whose reports agree with the target may round to zero", {
  value <- risk_diff(c(1, 1 + 2^-36), c(0, 0), 0, "MLE")$smse_boot
  expect_true(is.finite(value) && value >= 0 && value < 1e-100)
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

normalise_document_text <- function(text) {
  text <- paste(text, collapse = " ")
  text <- gsub(paste0("\\\\code[[:space:]]*\\{[[:space:]]*",
                      "((?:[^{}\\\\]|\\\\.)*?)[[:space:]]*\\}"),
               "\\1", text, perl = TRUE)
  text <- gsub("\\\\([%{}^])", "\\1", text, perl = TRUE)
  text <- gsub("`", "", text, fixed = TRUE)
  gsub("[[:space:]]+", " ", text)
}

# A sentence of a document, checked as a whole: the fixed text from its first
# word to its full stop, at the start of the text or after the full stop of the
# sentence before it. A sentence with words put in front of it, or with words
# between its last word and its full stop, does not match.
has_document_sentence <- function(text, sentence) {
  at <- gregexpr(sentence, text, fixed = TRUE)[[1]]
  if (at[1] < 0) return(FALSE)
  any(vapply(at, function(i) i == 1 || identical(substr(text, i - 2, i - 1), ". "),
             logical(1)))
}

# The sentences on the limits of the double range, as README.md and the help
# pages state them after markup is removed.
risk_limit_sentences <- c(
  terms = paste0(
    "The four scaled risks (smse_info, smse_boot, diff, smse_select) are summed on ",
    "the log scale: each term is formed as a whole from the logarithm of its weight ",
    "and of its two differences, so the sum never forms a weight, a squared ",
    "difference or a partial product on its own, and none of them can stop a call ",
    "or remove a term by leaving the double range."),
  checked = paste0(
    "Apart from requiring that the log weights and differences it receives are ",
    "finite, the sum checks the range only on the logarithm of the whole sum and ",
    "on the value rebuilt from it."),
  stop = paste0(
    "A risk whose logarithm, as this summation computes it, is below log(2^-1074) ",
    "stops the call with representable and the field name, even when ordinary ",
    "arithmetic would round the value up to 2^-1074; a value above the largest ",
    "double stops the call in the same way."),
  follow = paste0(
    "Both decisions follow the logarithm computed from the reports, the target and ",
    "the log weights in double precision, not the exact risk, and the two can lie ",
    "on different sides of 2^-1074 or of the largest double."),
  not_close = paste0(
    "This is not confined to risks close to those limits: when a report and the ",
    "target nearly agree, their difference keeps few significant digits and the ",
    "computed risk can be wrong by a large factor."),
  call_stopped = paste0(
    "In risk_diff(rep(1000, 3), c(-2 * r, 0, 2 * r), -r + 2^-52, \"MLE\") with ",
    "r = 0.685794632511088 the report and the target agree to a relative 1.1e-13, ",
    "the exact smse_boot is 2.05 times 2^-1074, which is representable, and the ",
    "call stopped when it was measured (R 4.6.1, aarch64 macOS)."),
  call_returned = paste0(
    "With r = 0.6907256300571373 and the ability -r + 2^-49 the exact value is ",
    "0.96 times 2^-1074, which is below the smallest positive double, and the call ",
    "returned 2^-1074."),
  no_claim = paste0(
    "Beyond the calls named in the specification, the package does not claim that ",
    "the outcome agrees with the side of either limit on which the exact risk lies."),
  diff_stop = paste0(
    "A nonzero diff whose logarithm of the absolute value, computed in the same ",
    "way, is below log(2^-1074) stops the call with representable and the field ",
    "name."))

# The attribution of the error of the order of 1e-13, through the words after
# the number.
risk_log_rounding_clause <- paste0(
  "Rounding that logarithm alone contributes a relative error of the order of ",
  "1e-13 where the logarithm is near 700 in absolute value;")

# Sentences of earlier versions that state more than the package does.
retired_risk_sentences <- c(
  "near the ends of the double range its relative error is of the order of 1e-13",
  "within a relative 1e-12 of 2^-1074",
  "inside the band at the two ends of the double range",
  "A positive risk below 2^-1074 stops the call",
  "A nonzero diff with absolute value below 2^-1074 stops the call",
  "A risk that is itself a subnormal number")

# Requirement: Numerical range
# Scenario: Documents state what a risk difference of 0 means
# spec-example: irtse-standard-errors / documents on the limits of the double range
test_that("README and risk help explain a zero diff as not distinguished", {
  texts <- list()
  readme <- testthat::test_path("..", "..", "README.md")
  if (file.exists(readme)) texts$README <- readLines(readme, warn = FALSE)
  rd_path <- testthat::test_path("..", "..", "man", "risk_diff.Rd")
  # In an installed package the page is read from the Rd database. Its pieces
  # are joined without a separator, which gives the text of the source file; a
  # space between the pieces would stand between a code span and the comma,
  # semicolon or bracket next to it.
  texts$risk_diff <- if (file.exists(rd_path)) readLines(rd_path, warn = FALSE) else {
    paste(as.character(tools::Rd_db("irtse")[["risk_diff.Rd"]]), collapse = "")
  }
  for (name in names(texts)) {
    text <- gsub("\\s+", " ", paste(texts[[name]], collapse = " "))
    expect_match(text, "diff[^.]*not distinguished[^.]*double precision",
                 info = name)
    text <- normalise_document_text(texts[[name]])
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
})

# Requirement: Numerical range
# Scenario: Documents state what a risk difference of 0 means
test_that("The document sentence check takes whole sentences and reads through markup", {
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
# Scenario: Documents state what a risk difference of 0 means
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
  path <- testthat::test_path("..", "..", ".Rbuildignore")
  if (!file.exists(path)) skip("Build ignore file not available in installed package")
  patterns <- readLines(path, warn = FALSE)
  expect_true(any(vapply(patterns, function(pattern) {
    grepl(pattern, ".git", perl = TRUE, ignore.case = TRUE)
  }, logical(1))))
})
