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

# Requirement: Estimators on the sufficient statistic
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
test_that("Shifted WLE designs", {
  shifted <- tryCatch(se_boot(1, 1000, 1000, "WLE"), error = identity)
  expect_type(shifted, "double")
  if (is.numeric(shifted)) expect_lt(abs(shifted / log(3) - 1), 1e-12)
  steep <- tryCatch(se_boot(rep(1000, 2), c(0, 0), 0, "WLE"), error = identity)
  expect_true(is.numeric(steep) && is.finite(steep) && steep > 0)
})

# Requirement: Estimators on the sufficient statistic
test_that("Two WLE roots inside one grid cell", {
  expect_error(
    se_boot(c(40, 41, 3.480672), c(0, 0.1, 0), 0, "WLE"),
    "T = 0", fixed = TRUE
  )
})

# Requirement: Numerical range
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
    list(call = function() delta2(1e-30, 0, 0, "MLE"), expected = 0,
         must_return = TRUE),
    list(call = function() delta3(1e-30, 0, 0, "MLE"),
         expected = -121 / (3 * 1e-60), tolerance = 1e-10, must_return = TRUE),
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
  forbidden_path <- paste0("/", "Users", "/")
  private_issue <- paste0("Psych", "Quant/[^\r\n]*", "#", "[0-9]+")
  bad <- vapply(unique(paths), function(path) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    length(grepRaw(forbidden_path, bytes, fixed = TRUE)) > 0L ||
      length(grepRaw(private_issue, bytes)) > 0L
  }, logical(1))
  expect_false(any(bad), info = paste(names(bad)[bad], collapse = ", "))
})

# Requirement: Numerical range
test_that("Log information and cumulants survive unrepresentable intermediates", {
  expect_lt(abs(se_info(rep(1, 6), rep(0, 6), 800) /
                  (exp(400) / sqrt(6)) - 1), 1e-10)
  # eta = 1000: a^r overflows for r >= 6, while a^r P Q is finite.
  values <- tryCatch(design_cumulants(1e60, 0, 1e-57), error = identity)
  if (inherits(values, "error")) {
    expect_match(conditionMessage(values), "representable", fixed = TRUE)
  } else {
    expected <- (-1)^(2:7) * exp((2:7) * log(1e60) - 1000)
    expect_true(all(is.finite(values)))
    expect_true(all(values != 0))
    expect_lt(max(abs(values / expected - 1)), 1e-10)
  }
})

# Requirement: Numerical range
test_that("A saturated item does not erase representable coefficients", {
  theta <- qlogis(1/5)
  # Only the first of two items contributes: k_r is halved, multiplying
  # Delta_2 by 2^3 and Delta_3 by 2^4 relative to the rational one-item case.
  expect_lt(abs(delta2(c(1, 1e100), c(0, 1), theta, "MLE") /
                  (-43425/64) - 1), 1e-10)
  expect_lt(abs(delta2(c(1, 1e100), c(0, 1), theta, "WLE") /
                  (-15075/64) - 1), 1e-10)
  expect_lt(abs(delta3(c(1, 1e100), c(0, 1), theta, "MLE") /
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
test_that("Underflowed score tails cannot certify a false MLE root", {
  expect_error(se_boot(c(1, 2), c(-1000, 1000), 3000, "MLE"),
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
test_that("Extreme range rows require their specified values and errors", {
  expect_error(se_boot(c(1, 2), c(-1000, 1000), 3000, "MLE"),
               "representable", fixed = TRUE)
  expect_lt(abs(delta3(1, 0, 75, "MLE") / (-(9689/768) * exp(300)) - 1),
            1e-10)
})

# Requirement: Numerical range; Information-based standard error
test_that("Information does not require representable higher cumulants", {
  expect_lt(abs(se_info(c(1, 1), c(0, 1000), 0) / 2 - 1), 1e-12)
})

# Requirement: Numerical range; Exact scaled risk of the selection rule
test_that("Risk evaluates both MLE coefficients even outside the margin", {
  expect_error(risk_diff(rep(1e-153, 3), rep(0, 3), 0, "MLE", margin = 0),
               "representable", fixed = TRUE)
  expect_error(se_select(rep(1e-153, 3), rep(0, 3), -log(2) / 1e-153,
                         "MLE", margin = 0), "representable", fixed = TRUE)
})

# Requirement: Numerical range
test_that("Small logits keep their odd cumulants", {
  expect_lt(abs(design_cumulants(1, 0, 1e-18)[["k3"]] / -1.25e-19 - 1), 1e-10)
  expect_lt(abs(delta2(1, 0, 1e-18, "MLE") / -1.4e-35 - 1), 1e-10)
  expect_lt(abs(delta2(1, 0, 1e-18, "WLE") / -5e-36 - 1), 1e-10)
  expect_lt(abs(delta2(1e-20, 0, 1, "MLE") / -14 - 1), 1e-10)
  expect_lt(abs(delta2(1e-20, 0, 1, "WLE") / -5 - 1), 1e-10)
  selected <- se_select(rep(1e-16, 2), rep(0, 2), 0.1, "WLE", margin = 0)
  expect_identical(selected$choice, "info")
  expect_lt(abs(selected$delta2 / -0.05 - 1), 1e-10)
  expect_lt(abs(se_diagnostic(rep(1, 6), rep(0, 6), 1e-18,
                             "WLE")$sign_delta2 / -1 - 1), 1e-10)
})

# Requirement: Numerical range
test_that("Large common shifts preserve law-based reports", {
  expect_lt(abs(se_boot(1, 1e16, 1e16, "WLE") / log(3) - 1), 1e-12)
  r <- risk_diff(1, 1e16, 1e16, "WLE", margin = 0)
  expect_lt(abs(r$h_true / log(3) - 1), 1e-10)
  expect_lt(abs(r$smse_info / (4 / sqrt(3) - log(3))^2 - 1), 1e-10)
  expect_lt(abs(r$smse_boot / ((sqrt(3) / 2 - 1) * log(3))^2 - 1), 1e-10)
  expect_lt(abs(r$smse_select / (4 / sqrt(3) - log(3))^2 - 1), 1e-10)
  expect_lt(abs(se_select(1, 1e16, 1e16, "WLE", margin = 0)$se_boot /
                  log(3) - 1), 1e-12)
})

# Requirement: Numerical range
test_that("Cancellation between items removes only exact opposite pairs", {
  expected <- plogis(-40) * plogis(40) * (plogis(40) - plogis(-40)) / 3
  for (b in list(c(-1, 40, 1), c(-1, 1, 40))) {
    expect_lt(abs(design_cumulants(rep(1, 3), b, 0)[["k3"]] / expected - 1),
              1e-12)
  }
  value <- tryCatch(design_cumulants(rep(2^100, 2),
                                    c(1, -(1 + 2^-50)) / 2^100, 0),
                    error = identity)
  if (inherits(value, "error")) {
    expect_match(conditionMessage(value), "not resolved", fixed = TRUE)
  } else {
    expect_true(is.finite(value[["k3"]]) && value[["k3"]] != 0)
  }
})

# Requirement: Input validation at the boundary
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
  patterns <- c(paste0("/", "Users", "/"), paste0("/", "home", "/"),
                paste0("[A-Za-z]", ":", "\\\\"),
                paste0("Psych", "Quant/[^\r\n]*", "#", "[0-9]+"),
                paste0("/", "issues", "/", "[0-9]+"))
  bad <- vapply(paths, function(path) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    any(vapply(patterns, function(pattern) length(grepRaw(pattern, bytes)) > 0L,
               logical(1)))
  }, logical(1))
  expect_false(any(bad), info = paste(paths[bad], collapse = ", "))
})

# Requirement: Exact law of the sufficient statistic
test_that("Nearly equal discriminations are not merged", {
  d <- 2^-40
  expected <- 2 * d / (1 + (1 + d)^2)
  expect_lt(abs(se_boot(c(1, 1 + d), c(0, 0), 0, "MLE") / expected - 1), 1e-6)
  expect_null(irtse:::lattice_weights(c(1, 1 + d)))
})

# Requirement: Exact law of the sufficient statistic; Limits of exact enumeration
test_that("Representation tolerance preserves decimal and earlier lattices", {
  expect_lt(abs(se_boot(c(0.1, 0.2, 0.3), c(0, 0, 0), 0, "WLE") /
                  (10 * se_boot(c(1, 2, 3), c(0, 0, 0), 0, "WLE")) - 1), 1e-10)
  cases <- list(c(0.1, 0.2, 0.3), c(1, 1000001/1000), rep(c(1, 997), 600),
                c(1, 999), c(1, 2), c(1, 2, 3, 1.5), rep(2, 5))
  for (a in cases) expect_false(is.null(irtse:::lattice_weights(a)))
})

# Requirement: Estimators on the sufficient statistic
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
