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

rule_coefficients <- function(a, eta, estimator) {
  design <- design_ratios(a, eta)
  list(delta2 = delta2_from_ratios(design, estimator),
       delta3 = if (estimator == "MLE") delta3_from_ratios(design) else NA_real_)
}

selection_at <- function(a, eta, estimator, margin) {
  coefficients <- rule_coefficients(a, eta, estimator)
  selection_decision(coefficients$delta2, coefficients$delta3, estimator, margin)$choice
}
