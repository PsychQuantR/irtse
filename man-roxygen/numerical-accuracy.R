#' @section Numerical accuracy:
#' Results depend on `theta` and `b` only through `a * (theta - b)` computed in
#' double precision; computations using the law of T and diagnostic zero searches
#' use coordinates relative to the midpoint of the difficulties.
#' A design average whose nonzero terms cancel to zero after rounding stops with
#' `not resolved`; accuracy degrades with cancellation between items and near
#' the zeros of a coefficient.
#' Log scales and normalised cumulant ratios avoid intermediate overflow and
#' underflow, while a nonzero result outside double precision stops with
#' `representable` (except underflowed probabilities in `risk_diff()`'s `p_select`,
#' which are returned as zero).
