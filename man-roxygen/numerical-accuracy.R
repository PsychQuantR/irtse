#' @section Supported range:
#' Every exported function requires `1e-3 <= a_j <= 1e3`, `abs(b_j) <= 1e3`
#' and `abs(theta) <= 1e3`, inclusive. The bounds are the doubles represented
#' by the R literals `1e-3` and `1e3`. This range binds the supplied arguments only.
#' Estimates and internal evaluation points may lie outside it and are neither
#' rejected nor clamped. This includes root intervals and the window, grid points
#' and `nearest_zero` of `se_diagnostic()`.
#' With `margin`, `risk_diff()` evaluates the selection rule at every estimate,
#' including estimates outside the range; supplying such an estimate as `theta`
#' to `se_select()` is rejected.
#' Numeric arguments (`a`, `b`, `theta`, `margin`) that carry a class are rejected
#' before arithmetic; plain numeric vectors with names are accepted and their
#' names are ignored.
#'
#' @section Numerical accuracy:
#' At the supplied ability, information and coefficients use
#' `a * (theta - b)` computed from the supplied inputs in double precision.
#' The law of T uses these logits with the identified discriminations described
#' below. In law-based computations,
#' estimates remain offsets from the midpoint of the difficulties, while
#' diagnostic zero searches use offsets from the supplied `theta`.
#' Accuracy depends on cancellation between items, each item's cumulant polynomial
#' near its zeros, the subtraction `theta - b` of nearly equal numbers,
#' and proximity to a coefficient's
#' zeros; a design average whose nonzero terms cancel to zero after rounding stops
#' with `not resolved`, except in the zero search of `se_diagnostic()`, where an
#' unresolved average counts as a zero.
#' Log scales and normalised cumulant ratios reduce intermediate overflow and
#' underflow. A call can still stop with `representable` when an intermediate
#' value, such as a scaled cumulant ratio or a squared report before its weight
#' is applied, is not representable, even if the returned quantity would be.
#' Weighted score terms are formed on the log scale, including the discrimination.
#' The function and its derivative use the same scale at each point for signs,
#' exact-zero tests and Newton steps, so the location of a root does not depend
#' on the double range, even when every score term is smaller than the smallest
#' positive double.
#' Inputs within the supported range can still yield a nonzero result
#' outside double precision, as in `delta2(1, 0, 300, "MLE")`, or a nonzero logit
#' whose `tanh(eta / 2)` is zero, as in `delta2(1, 0, 2^-1074, "MLE")`.
#' These cases, and any non-finite logit, stop with `representable`
#' (underflowed probabilities in `risk_diff()`'s `p_select` are returned as zero).
#' Discriminations whose ratios satisfy the lattice tolerance of eight relative
#' machine epsilons are identified with the commensurable design, even when
#' the supplied doubles are distinct. The identified design is used throughout
#' the law-based functions: totals, estimating equations, estimates, the law at
#' `theta` and at each estimate, and every report and rule coefficient at an
#' estimate. The `se_info`, `delta2` and `delta3` fields of `se_select()` still use
#' the supplied discriminations. Identified totals form one atom and can give
#' an exactly zero bootstrap SE.
#' A third exception to nonzero results is a risk term whose report agrees with
#' its target to within rounding: its squared difference is zero. For example,
#' `risk_diff(c(1, 1 + 2^-36), c(0, 0), 0, "MLE")$smse_boot` can be zero even
#' though the mathematical risk is positive. A positive risk that underflows
#' stops the call with `representable`. A risk that is itself a subnormal number
#' has few significant bits. For `diff`, a term is zero when either factor in
#' `w * (A - B) * (A + B - 2 * h)` is exactly zero; cancelling positive and negative
#' parts can also give zero, but underflow of nonzero terms cannot.
