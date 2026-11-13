validate_design <- function(a, b, theta) {
  if (!is.numeric(a) || !is.null(dim(a)) || !length(a) ||
      any(!is.finite(a)) || any(a <= 0)) {
    stop("a must be a non-empty numeric vector of finite positive values.",
         call. = FALSE)
  }
  if (!is.numeric(b) || !is.null(dim(b)) || length(b) != length(a) ||
      any(!is.finite(b))) {
    stop("b must be a finite numeric vector with the same length as a.",
         call. = FALSE)
  }
  if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta)) {
    stop("theta must be a single finite number.", call. = FALSE)
  }
  invisible(NULL)
}

validate_estimator <- function(estimator) {
  if (missing(estimator) || !is.character(estimator) ||
      length(estimator) != 1L || is.na(estimator) ||
      !estimator %in% c("MLE", "WLE")) {
    stop('estimator must be supplied as exactly "MLE" or "WLE".',
         call. = FALSE)
  }
  invisible(NULL)
}

validate_margin <- function(margin) {
  if (missing(margin) || !is.numeric(margin) || length(margin) != 1L ||
      !is.finite(margin) || margin < 0) {
    stop("margin must be supplied as a single finite non-negative number.",
         call. = FALSE)
  }
  invisible(NULL)
}

validate_law_estimator <- function(a, estimator) {
  validate_estimator(estimator)
  if (estimator == "MLE" && length(a) < 2L) {
    stop("MLE requires at least two items for non-extreme responses.",
         call. = FALSE)
  }
  invisible(NULL)
}
