#' @section Numerical accuracy:
#' At the supplied ability, information, coefficients and the law of T use
#' `a * (theta - b)` computed from the supplied inputs in double precision;
#' estimates remain offsets from the midpoint of the difficulties, while
#' diagnostic zero searches use offsets from the supplied `theta`.
#' Accuracy depends on cancellation between items, each item's cumulant polynomial
#' near its zeros, the subtraction `theta - b`, and proximity to a coefficient's
#' zeros; a design average whose nonzero terms cancel to zero after rounding stops
#' with `not resolved`, except in the zero search of `se_diagnostic()`, where an
#' unresolved average counts as a zero.
#' Log scales and normalised cumulant ratios avoid intermediate overflow and
#' underflow, while a nonzero result outside double precision, a non-finite logit,
#' or a nonzero logit whose `tanh(eta / 2)` is zero stops with `representable`
#' (underflowed probabilities in `risk_diff()`'s `p_select` are returned as zero).
#' Discriminations within 8 units in the last place of a commensurable design are
#' treated as that design, so totals within this lattice tolerance are identified
#' as one atom and can give an exactly zero bootstrap SE.
