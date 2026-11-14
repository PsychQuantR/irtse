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

# Requirement: Input validation at the boundary
test_that("Numeric arguments carrying a class are rejected at the boundary", {
  error_message <- function(expr) tryCatch({ force(expr); "" }, error = conditionMessage)
  cases <- list(
    list(call = function() se_info(stats::ts(c(1, 1), frequency = 2),
                                   stats::ts(c(0, 0)), 0), prefix = "^a must"),
    list(call = function() se_info(c(1, 1), stats::ts(c(0, 0)), 0),
         prefix = "^b must"),
    list(call = function() se_info(c(1, 1), c(0, 0),
                                   structure(0, class = "units")), prefix = "^theta must"),
    list(call = function() se_select(c(1, 1), c(0, 0), 0, "WLE",
                                     margin = structure(0, class = "units")),
         prefix = "^margin must"),
    list(call = function() se_boot(structure(c(1, 1), class = "integer64"),
                                   c(0, 0), 0, "WLE"), prefix = "^a must")
  )
  for (case in cases) expect_match(error_message(case$call()), case$prefix)
  expect_lt(abs(se_info(c(x = 1, y = 1), c(u = 0, v = 0), c(t = 0)) /
                  sqrt(2) - 1), 1e-12)
})

# Requirement: Input validation at the boundary
# Scenario: Numeric arguments that carry a class
test_that("A classed estimator is rejected even without class methods", {
  expect_error(se_boot(c(1, 1), c(0, 0), 0,
                       structure("WLE", class = "irtse_plain")), "^estimator must")
})

with_trap_methods <- function(code) {
  generics <- c("is.numeric", "is.character", "length", "dim", "is.finite",
                "is.na", "as.character", "Ops", "Math", "Summary")
  if (exists("mtfrm", baseenv())) generics <- c(generics, "mtfrm")
  method_names <- paste0(generics, ".irtse_trap")
  present <- vapply(method_names, exists, logical(1),
                    envir = .GlobalEnv, inherits = FALSE)
  previous <- mget(method_names[present], envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    rm(list = method_names, envir = .GlobalEnv)
    list2env(previous, envir = .GlobalEnv)
  })
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  for (method in method_names) {
    assign(method, function(...) {
      counter$calls <- counter$calls + 1L
      stop("irtse trap method invoked", call. = FALSE)
    }, envir = .GlobalEnv)
  }
  code(counter)
}

# Requirement: Input validation at the boundary
# Scenario: A class whose methods raise errors
test_that("All five arguments reject a live trap before any method runs", {
  with_trap_methods(function(counter) {
    x <- structure(c(1, 1), class = "irtse_trap")
    message <- tryCatch(is.numeric(x), error = conditionMessage)
    expect_identical(message, "irtse trap method invoked")
    expect_identical(counter$calls, 1L)
    plain <- list(a = c(1, 1), b = c(0, 0), theta = 0, estimator = "WLE", margin = 0)
    for (argument in names(plain)) {
      args <- plain
      args[[argument]] <- structure(if (argument == "b") c(1, 1) else args[[argument]],
                                    class = "irtse_trap")
      counter$calls <- 0L
      message <- tryCatch({ do.call(se_select, args); "" }, error = conditionMessage)
      expect_match(message, paste0("^", argument, " must"))
      expect_identical(counter$calls, 0L, info = argument)
    }
  })
})

# Requirement: Input validation at the boundary
# Scenario: A class whose methods raise errors
test_that("Every exported function rejects a trapped a before dispatch", {
  with_trap_methods(function(counter) {
    x <- structure(c(1, 1), class = "irtse_trap")
    message <- tryCatch(is.numeric(x), error = conditionMessage)
    expect_identical(message, "irtse trap method invoked")
    expect_identical(counter$calls, 1L)
    for (name in c("se_info", "se_boot", "design_cumulants", "delta2", "delta3",
                   "risk_diff", "se_diagnostic", "se_select")) {
      fun <- getExportedValue("irtse", name)
      args <- list(a = x, b = c(0, 0), theta = 0)
      if ("estimator" %in% names(formals(fun))) args$estimator <- "WLE"
      if ("margin" %in% names(formals(fun))) args$margin <- 0
      counter$calls <- 0L
      message <- tryCatch({ do.call(fun, args); "" }, error = conditionMessage)
      expect_match(message, "^a must")
      expect_identical(counter$calls, 0L, info = name)
    }
  })
})

# Requirement: Input validation at the boundary
test_that("Every exported function rejects classed numeric arguments in order", {
  error_message <- function(expr) tryCatch({ force(expr); "" }, error = conditionMessage)
  plain <- list(a = c(1, 1), b = c(0, 0), theta = 0)
  for (name in c("se_info", "se_boot", "design_cumulants", "delta2", "delta3",
                 "risk_diff", "se_diagnostic", "se_select")) {
    fun <- getExportedValue("irtse", name)
    args <- plain
    if ("estimator" %in% names(formals(fun))) args$estimator <- "MLE"
    if ("margin" %in% names(formals(fun))) args$margin <- 0
    for (argument in intersect(c("a", "b", "theta", "margin"), names(args))) {
      invalid <- args
      invalid[[argument]] <- structure(invalid[[argument]], class = "units")
      expect_match(error_message(do.call(fun, invalid)),
                   paste0("^", argument, " must"), info = paste(name, argument))
    }
    args$a <- stats::ts(args$a, frequency = 2)
    args$b <- stats::ts(args$b)
    expect_match(error_message(do.call(fun, args)), "^a must", info = name)
  }
})
