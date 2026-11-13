# Requirement: Input validation at the boundary
test_that("Supported range accepts each inclusive bound", {
  cases <- list(
    list(args = list(1e-3, 0, 0), expected = 2000),
    list(args = list(1e3, 0, 0), expected = 0.002),
    list(args = list(1, 1e3, 0), expected = exp(500) * (1 + exp(-1000))),
    list(args = list(1, -1e3, 0), expected = exp(500) * (1 + exp(-1000))),
    list(args = list(1, 0, 1e3), expected = exp(500) * (1 + exp(-1000))),
    list(args = list(1, 0, -1e3), expected = exp(500) * (1 + exp(-1000)))
  )
  for (case in cases) {
    expect_lt(abs(do.call(se_info, case$args) / case$expected - 1), 1e-12)
  }
})

# Requirement: Input validation at the boundary
test_that("Supported range rejects the adjacent doubles and states the bounds", {
  cases <- list(
    list(args = list(1e-3 - 2^-62, 0, 0), argument = "a"),
    list(args = list(1e3 + 2^-43, 0, 0), argument = "a"),
    list(args = list(1, 1e3 + 2^-43, 0), argument = "b"),
    list(args = list(1, -(1e3 + 2^-43), 0), argument = "b"),
    list(args = list(1, 0, 1e3 + 2^-43), argument = "theta"),
    list(args = list(1, 0, -(1e3 + 2^-43)), argument = "theta")
  )
  for (case in cases) {
    message <- tryCatch({
      do.call(se_info, case$args)
      ""
    }, error = conditionMessage)
    expect_match(message, paste0("^", case$argument, " must"))
    if (case$argument == "a") expect_match(message, "1e-3", fixed = TRUE)
    expect_match(message, "1e3", fixed = TRUE)
  }
})

# Requirement: Numerical range
test_that("Inputs outside the supported range name the supplied argument", {
  cases <- list(
    list(call = function() delta2(1e-20, 0, 1, "MLE"), prefix = "^a must"),
    list(call = function() se_boot(1, 1e16, 1e16, "WLE"), prefix = "^b must"),
    list(call = function() delta2(c(1, 1e308), c(0, 10), 0, "MLE"),
         prefix = "^a must"),
    list(call = function() delta2(2^-1074, 0, 1, "MLE"), prefix = "^a must"),
    list(call = function() se_boot(c(1, 2), c(-1000, 1000), 3000, "MLE"),
         prefix = "^theta must")
  )
  for (case in cases) {
    message <- tryCatch({ case$call(); "" }, error = conditionMessage)
    expect_match(message, case$prefix)
  }
})

# Requirement: Numerical range; Exact scaled risk of the selection rule
test_that("Estimates outside the supported range remain valid internal points", {
  expect_lt(abs(se_boot(1e-3, 0, 0, "WLE") / (1000 * log(3)) - 1), 1e-12)
  risk <- risk_diff(1e-3, 0, 0, "WLE", margin = 0)
  expect_lt(abs(risk$h_true / (1000 * log(3)) - 1), 1e-10)
  expect_lt(abs(risk$smse_info / (1e6 * (4 / sqrt(3) - log(3))^2) - 1),
            1e-10)
  expect_lt(abs(risk$smse_boot /
                  (1e6 * ((sqrt(3) / 2 - 1) * log(3))^2) - 1), 1e-10)
  expect_lt(abs(risk$smse_select / (1e6 * (4 / sqrt(3) - log(3))^2) - 1),
            1e-10)
  expect_identical(risk$smse_select, risk$smse_info)
  expect_equal(risk$p_select, c(info = 1, boot = 0, none = 0), tolerance = 1e-12)
  expect_error(se_select(1e-3, 0, 1000 * log(3), "WLE", margin = 0),
               "^theta must")
})

# Requirement: Input validation at the boundary
test_that("Supported range validation covers every exported function", {
  calls <- list(
    se_info = function() se_info(c(1, 2000), c(0, 0), 0),
    se_boot = function() se_boot(c(1, 2000), c(0, 0), 0, "MLE"),
    design_cumulants = function() design_cumulants(c(1, 2000), c(0, 0), 0),
    delta2 = function() delta2(c(1, 2000), c(0, 0), 0, "MLE"),
    delta3 = function() delta3(c(1, 2000), c(0, 0), 0, "MLE"),
    risk_diff = function() risk_diff(c(1, 2000), c(0, 0), 0, "MLE", margin = 0),
    se_diagnostic = function() se_diagnostic(c(1, 2000), c(0, 0), 0, "MLE"),
    se_select = function() se_select(c(1, 2000), c(0, 0), 0, "MLE", margin = 0)
  )
  for (name in names(calls)) {
    message <- tryCatch({ calls[[name]](); "" }, error = conditionMessage)
    expect_match(message, "^a must", info = name)
  }
})

# Requirement: Input validation at the boundary
test_that("Range validation follows a then b then theta before other arguments", {
  cases <- list(
    list(call = function() se_info(c(1, 2000), c(0, 2000), 0),
         prefix = "^a must"),
    list(call = function() se_info(c(1, 1), c(0, 2000), 2000),
         prefix = "^b must"),
    list(call = function() se_boot(c(1, 1), c(0, 0), 2000),
         prefix = "^theta must"),
    list(call = function() se_select(c(1, 1), c(0, 0), 2000, "MLE", margin = -1),
         prefix = "^theta must")
  )
  for (case in cases) {
    message <- tryCatch({ case$call(); "" }, error = conditionMessage)
    expect_match(message, case$prefix)
  }
})

# Requirement: Package metadata and documentation
test_that("README states the supported range", {
  readme_path <- testthat::test_path("..", "..", "README.md")
  if (!file.exists(readme_path)) {
    skip("README not available in installed package")
  }
  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
  expect_match(readme, "1e-3", fixed = TRUE)
  expect_match(readme, "1e3", fixed = TRUE)
})

# Requirement: Package metadata and documentation
test_that("All eight help pages state the supported range and its scope", {
  for (name in c("se_info", "se_boot", "design_cumulants", "delta2", "delta3",
                 "risk_diff", "se_diagnostic", "se_select")) {
    rd_name <- paste0(name, ".Rd")
    rd_path <- testthat::test_path("..", "..", "man", rd_name)
    rd <- if (file.exists(rd_path)) readLines(rd_path, warn = FALSE) else {
      as.character(tools::Rd_db("irtse")[[rd_name]])
    }
    help_text <- gsub("\\s+", " ", paste(rd, collapse = " "))
    for (phrase in c("1e-3 <= a_j <= 1e3", "abs(b_j) <= 1e3",
                     "abs(theta) <= 1e3", "inclusive", "supplied arguments only",
                     "Estimates and internal evaluation points may lie outside")) {
      expect_match(help_text, phrase, fixed = TRUE, info = name)
    }
  }
})
