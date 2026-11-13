selection_decision <- function(d2, d3, estimator, margin) {
  finite_value(d2, "Delta_2")
  within <- abs(d2) <= margin
  if (within && estimator == "WLE") {
    return(list(choice = "none", within_margin = TRUE,
                reason = paste("No choice is made within the Delta_2 margin",
                               "because WLE has no Delta_3 coefficient.")))
  }
  coefficient <- if (within) d3 else d2
  finite_value(coefficient, "selection coefficient")
  label <- if (within) "Delta_3" else "Delta_2"
  if (coefficient == 0) {
    return(list(choice = "none", within_margin = within,
                reason = paste("No choice is made because", label, "is zero.")))
  }
  negative <- coefficient < 0
  choice <- if (negative) "info" else "boot"
  direction <- if (negative) "negative" else "positive"
  list(choice = choice, within_margin = within,
       reason = sprintf("%s is %s, so the rule chooses the %s SE.",
                        label, direction,
                        if (negative) "information" else "bootstrap"))
}

selection_at <- function(a, b, theta, estimator, margin) {
  d2 <- delta2(a, b, theta, estimator)
  within <- abs(d2) <= margin
  d3 <- if (within && estimator == "MLE") {
    delta3(a, b, theta, estimator)
  } else {
    NA_real_
  }
  selection_decision(d2, d3, estimator, margin)$choice
}
