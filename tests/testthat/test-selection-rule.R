# Requirement: Selection by the second-order coefficient outside the margin
test_that("Decisions outside the margin", {
  cases <- list(
    list(a = rep(1, 6), b = rep(0, 6), theta = qlogis(3/10),
         estimator = "MLE", expected_delta2 = -48800/3087, choice = "info"),
    list(a = c(1, 1, 1), b = c(-log(4), -log(4), log(4)), theta = 0,
         estimator = "MLE", expected_delta2 = 2175/512, choice = "boot"),
    list(a = rep(1, 6), b = rep(0, 6), theta = qlogis(1/5),
         estimator = "WLE", expected_delta2 = -15075/512, choice = "info")
  )

  for (case in cases) {
    selection <- se_select(
      a = case$a, b = case$b, theta = case$theta, estimator = case$estimator,
      margin = 1
    )
    info <- se_info(a = case$a, b = case$b, theta = case$theta)
    boot <- se_boot(
      a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
    )

    expect_type(selection, "list")
    expect_named(
      selection,
      c("choice", "se", "se_info", "se_boot", "delta2", "delta3",
        "within_margin", "reason"),
      ignore.order = TRUE
    )
    expect_identical(selection$choice, case$choice)
    expect_false(selection$within_margin)
    expect_identical(selection$se, if (case$choice == "info") info else boot)
    expect_identical(selection$se_info, info)
    expect_identical(selection$se_boot, boot)
    expect_identical(
      selection$delta2,
      delta2(
        a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
      )
    )
    expect_equal(selection$delta2, case$expected_delta2, tolerance = 1e-10)
    if (case$estimator == "MLE") {
      expect_identical(
        selection$delta3,
        delta3(a = case$a, b = case$b, theta = case$theta, estimator = "MLE")
      )
    } else {
      expect_true(is.na(selection$delta3))
      expect_type(selection$delta3, "double")
    }
    expect_match(
      selection$reason,
      "(?i)(delta[ _-]?2|Δ[ _-]?[2₂]|second[ -]order|二階)", perl = TRUE
    )
  }
})

# Requirement: Third-order coefficient inside the margin
test_that("Decisions inside the margin", {
  cases <- list(
    list(a = rep(1, 6), b = rep(0, 6), theta = 0, estimator = "MLE",
         expected_delta3 = -121/3, choice = "info"),
    list(a = rep(1, 6), b = rep(c(log(5), -log(5)), 3), theta = 0,
         estimator = "MLE", expected_delta3 = 110889/625, choice = "boot"),
    list(a = rep(1, 6), b = rep(0, 6), theta = 0, estimator = "WLE",
         expected_delta3 = NA_real_, choice = "none")
  )

  for (case in cases) {
    selection <- se_select(
      a = case$a, b = case$b, theta = case$theta, estimator = case$estimator,
      margin = 1e-8
    )
    info <- se_info(a = case$a, b = case$b, theta = case$theta)
    boot <- se_boot(
      a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
    )

    expect_type(selection, "list")
    expect_named(
      selection,
      c("choice", "se", "se_info", "se_boot", "delta2", "delta3",
        "within_margin", "reason"),
      ignore.order = TRUE
    )
    expect_identical(selection$choice, case$choice)
    expect_true(selection$within_margin)
    if (case$estimator == "MLE") {
      expect_equal(selection$delta3, case$expected_delta3, tolerance = 1e-10)
      expect_identical(
        selection$delta3,
        delta3(a = case$a, b = case$b, theta = case$theta, estimator = "MLE")
      )
      expect_identical(selection$se, if (case$choice == "info") info else boot)
      expect_match(
        selection$reason,
        "(?i)(delta[ _-]?3|Δ[ _-]?[3₃]|third[ -]order|三階)", perl = TRUE
      )
    } else {
      expect_true(is.na(selection$delta3))
      expect_type(selection$delta3, "double")
      expect_true(is.na(selection$se))
      expect_type(selection$se, "double")
      expect_match(
        selection$reason,
        "(?i)(WLE|delta[ _-]?3|Δ[ _-]?[3₃]|third[ -]order|三階)", perl = TRUE
      )
    }
    expect_identical(selection$se_info, info)
    expect_identical(selection$se_boot, boot)
    expect_identical(
      selection$delta2,
      delta2(
        a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
      )
    )
  }
})

# Requirement: Explicit margin without a theoretical default
test_that("No choice for the WLE inside the margin and for a zero Delta_3", {
  r <- risk_diff(rep(1, 6), rep(0, 6), 0, "WLE", margin = 1e-8)

  expect_true(is.na(r$smse_select))
  expect_identical(names(r$p_select), c("info", "boot", "none"))
  expect_lt(abs(r$p_select[["info"]] - 0.6875), 1e-12)
  expect_lt(abs(r$p_select[["boot"]] - 0), 1e-12)
  expect_lt(abs(r$p_select[["none"]] - 0.3125), 1e-12)

  decision <- irtse:::selection_decision(0, 0, "MLE", 0)
  expect_identical(decision$choice, "none")
  expect_match(decision$reason, "Delta_3")
})

# Requirement: Explicit margin without a theoretical default
test_that("Invalid margins", {
  expect_error(se_select(rep(1, 6), rep(0, 6), 0, "MLE"), "margin")
  expect_error(
    se_select(rep(1, 6), rep(0, 6), 0, "MLE", margin = -1), "margin"
  )
  expect_error(
    se_select(rep(1, 6), rep(0, 6), 0, "MLE", margin = NA), "margin"
  )
  expect_error(
    se_select(rep(1, 6), rep(0, 6), 0, "MLE", margin = Inf), "margin"
  )
  expect_error(
    se_select(rep(1, 6), rep(0, 6), 0, "MLE", margin = c(1, 2)), "margin"
  )
  expect_error(
    risk_diff(rep(1, 6), rep(0, 6), 0, "MLE", margin = -1), "^margin must"
  )
})

# Requirement: Explicit margin without a theoretical default
test_that("Documented lack of guarantees", {
  getExportedValue("irtse", "se_select")
  rd_path <- testthat::test_path("..", "..", "man", "se_select.Rd")
  if (file.exists(rd_path)) {
    help_text <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")
  } else {
    help_text <- paste(
      as.character(tools::Rd_db("irtse")[["se_select.Rd"]]), collapse = ""
    )
  }

  expect_match(help_text, "provisional", fixed = TRUE)
  expect_match(help_text, "no guarantee", fixed = TRUE)
})

# Requirement: Documented scope of the selection rule
test_that("Help page states the scope", {
  getExportedValue("irtse", "se_diagnostic")
  rd_path <- testthat::test_path("..", "..", "man", "se_select.Rd")
  if (file.exists(rd_path)) {
    help_text <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")
  } else {
    help_text <- paste(
      as.character(tools::Rd_db("irtse")[["se_select.Rd"]]), collapse = ""
    )
  }

  for (phrase in c(
    "not a recommendation", "J^(-3/2)", "fixed test length", "not been refereed"
  )) {
    expect_match(help_text, phrase, fixed = TRUE)
  }
  expect_no_match(help_text, "recommended", fixed = TRUE)
})

# Requirement: Documented scope of the selection rule
test_that("README states the scope", {
  readme_path <- testthat::test_path("..", "..", "README.md")
  if (!file.exists(readme_path)) {
    skip("README not available in installed package")
  }

  getExportedValue("irtse", "se_diagnostic")
  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")

  for (phrase in c(
    "not a recommendation", "J^(-3/2)", "fixed test length", "not been refereed"
  )) {
    expect_match(readme, phrase, fixed = TRUE)
  }
  expect_no_match(readme, "recommended", fixed = TRUE)
})

# Requirement: Diagnostic without choosing
test_that("Diagnostic values", {
  cases <- list(
    list(a = rep(1, 6), b = rep(0, 6), theta = qlogis(3/10),
         estimator = "MLE", expected_delta2 = -48800/3087,
         expected_delta3 = -369968225/583443, expected_sign = -1,
         nearest_zero = 0),
    list(a = rep(1, 6), b = rep(0, 6), theta = qlogis(3/10),
         estimator = "WLE", expected_sign = -1, nearest_zero = 0),
    list(a = rep(c(1, 2), 25), b = rep(c(-1/2, 1), 25), theta = 0.03,
         estimator = "MLE", nearest_zero = 0.018053348)
  )

  for (case in cases) {
    if (case$estimator == "WLE") {
      expect_silent(diagnostic <- se_diagnostic(
        a = case$a, b = case$b, theta = case$theta, estimator = "WLE"
      ))
    } else {
      diagnostic <- se_diagnostic(
        a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
      )
    }

    expect_type(diagnostic, "list")
    expect_named(
      diagnostic,
      c("delta2", "delta3", "sign_delta2", "se_info", "nearest_zero",
        "zero_distance_se"),
      ignore.order = TRUE
    )
    if (is.null(case$expected_delta2)) {
      expect_identical(
        diagnostic$delta2,
        delta2(
          a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
        )
      )
    } else {
      expect_equal(diagnostic$delta2, case$expected_delta2, tolerance = 1e-10)
    }
    if (case$estimator == "WLE") {
      expect_true(is.na(diagnostic$delta3))
      expect_type(diagnostic$delta3, "double")
    } else if (is.null(case$expected_delta3)) {
      expect_identical(
        diagnostic$delta3,
        delta3(
          a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
        )
      )
    } else {
      expect_equal(diagnostic$delta3, case$expected_delta3, tolerance = 1e-10)
    }
    if (is.null(case$expected_sign)) {
      expect_equal(
        diagnostic$sign_delta2,
        sign(delta2(
          a = case$a, b = case$b, theta = case$theta, estimator = case$estimator
        ))
      )
    } else {
      expect_equal(diagnostic$sign_delta2, case$expected_sign)
    }
    expect_identical(
      diagnostic$se_info, se_info(a = case$a, b = case$b, theta = case$theta)
    )
    expect_lt(abs(diagnostic$nearest_zero - case$nearest_zero), 1e-8)
    expect_equal(
      diagnostic$zero_distance_se,
      abs(diagnostic$nearest_zero - case$theta) / diagnostic$se_info,
      tolerance = 1e-12
    )
  }
})

# Requirement: Diagnostic without choosing
test_that("No zero in the window", {
  diagnostic <- se_diagnostic(
    a = rep(1, 50), b = rep(0, 50), theta = qlogis(1/10), estimator = "MLE"
  )

  expect_true(is.na(diagnostic$nearest_zero))
  expect_type(diagnostic$nearest_zero, "double")
  expect_true(is.na(diagnostic$zero_distance_se))
  expect_type(diagnostic$zero_distance_se, "double")
})

# Requirement: Numerical range
test_that("Selection preserves the logits at the supplied ability", {
  a <- c(1, 1)
  b <- c(0, 80)
  theta <- 16 * exp(-80)
  selected <- se_select(a, b, theta, "MLE", margin = 0)
  expect_identical(selected$delta2, delta2(a, b, theta, "MLE"))
  expect_identical(selected$delta3, delta3(a, b, theta, "MLE"))
  expect_identical(selected$se_info, se_info(a, b, theta))
  expect_lt(abs(selected$delta2 / (5120 * exp(-160)) - 1), 1e-6)
  expect_identical(selected$choice, "boot")
})

# Requirement: Numerical range
test_that("Selection and standalone reports agree exactly on an ordinary design", {
  a <- c(1, 2, 1.5)
  b <- c(-0.5, 0.3, 1.2)
  theta <- 0.4
  for (estimator in c("MLE", "WLE")) {
    selected <- se_select(a, b, theta, estimator, margin = 0)
    expect_identical(selected$se_info, se_info(a, b, theta))
    expect_identical(selected$delta2, delta2(a, b, theta, estimator))
    expected_d3 <- if (estimator == "MLE") delta3(a, b, theta, estimator) else {
      suppressMessages(delta3(a, b, theta, estimator))
    }
    expect_identical(selected$delta3, expected_d3)
  }
})

# Requirement: Numerical range; Diagnostic without choosing
test_that("Diagnostic zero distances retain the supplied ability", {
  theta <- 0.1
  r <- se_diagnostic(c(1, 1e-3), c(0, 1e3), theta, "WLE")
  expect_true(is.finite(r$nearest_zero))
  expect_lt(abs(r$zero_distance_se -
                  abs(r$nearest_zero - theta) / r$se_info), 1e-9)
})

# Requirement: Numerical range; Diagnostic without choosing
test_that("Diagnostic offset search retains the known Rasch zero and distance", {
  theta <- qlogis(3/10)
  r <- se_diagnostic(rep(1, 6), rep(0, 6), theta, "MLE")
  expect_lt(abs(r$nearest_zero - 0), 1e-12)
  expect_lt(abs(r$zero_distance_se / (abs(qlogis(3/10)) / r$se_info) - 1),
            1e-12)
})
