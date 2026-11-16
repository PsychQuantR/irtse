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
#' Numeric arguments (`a`, `b`, `theta`, `margin`) and `estimator` that carry a
#' class are rejected before any class method is called; plain numeric vectors
#' with names are accepted and their names are ignored.
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
#' value, such as a scaled cumulant ratio, is not representable, even if the
#' returned quantity would be. The four scaled risks (`smse_info`, `smse_boot`,
#' `diff`, `smse_select`) are summed on the log scale: each term is formed as a
#' whole from the logarithm of its weight and of its two differences, so the sum
#' never forms a weight, a squared difference or a partial product on its own,
#' and none of them can stop a call or remove a term by leaving the double
#' range. Apart from requiring that the log weights and differences it receives
#' are finite, the sum checks the range only on the logarithm of the whole sum
#' and on the value rebuilt from it.
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
#' though the mathematical risk is positive. A risk whose logarithm, as this
#' summation computes it, is below `log(2^-1074)` stops the call with
#' `representable` and the field name, even when ordinary arithmetic would round
#' the value up to `2^-1074`; a risk whose value rebuilt from that logarithm is
#' not finite, or is 0, stops the call in the same way. The value of a scaled
#' risk is
#' rebuilt from its logarithm. Rounding that logarithm alone contributes a relative
#' error of the order of 1e-13 where the logarithm is near 700 in absolute value;
#' errors already in the reports and log weights add to it, and so does rounding
#' the result to a double. A returned risk in the subnormal range is a multiple
#' of `2^-1074`, so rounding the result to a subnormal double can give a much
#' larger
#' relative error: in `risk_diff(rep(1000, 3), c(-1.44, 0, 1.44), -0.72 + 5e-9,
#' "MLE")` the exact `smse_boot` is 1.6 times `2^-1074` and the call returned 2
#' times `2^-1074`. A
#' subnormal risk
#' is not claimed to be within `2^-1074` of the exact value, and no general
#' bound on the relative error of a returned risk is claimed.
#' Both decisions follow the logarithm computed from the reports, the target and
#' the log weights in double precision, not the exact risk, and the two can lie
#' on different sides of `2^-1074` or of the largest double. At `2^-1074` this
#' is not confined to risks close to the limit: when a report and the target
#' nearly
#' agree, their difference keeps few significant digits and the computed risk
#' can be wrong by a large factor. In `risk_diff(rep(1000, 3), c(-2 * r, 0, 2 *
#' r), -r + 2^-52, "MLE")` with `r` set to `0.685794632511088` the report and
#' the
#' target agree to a relative 1.1e-13, the exact `smse_boot` is 2.05 times
#' `2^-1074`, which is representable, and the call stopped. With `r` set to
#' `0.6907256300571373` and the ability `-r +
#' 2^-49` the exact value is 0.96 times `2^-1074`, which is below the smallest
#' positive double, and the call returned `2^-1074`. At the largest double, in
#' `risk_diff(rep(512, 2), c(-B, B), 0, "MLE")` with `B` set to `1042 * log(2) /
#' 512` the exact `smse_info` is above the largest double by a relative 9.9e-15
#' and the call returned 1.7976931348622732e308. In general the package does
#' not claim that the outcome agrees with the side of either limit on which the
#' exact risk lies. The outcomes of the calls named here were measured with R
#' 4.6.1 on aarch64 macOS and can differ on another platform.
#' For `diff`, cancellation multiplies the errors of its terms by its condition
#' number (the sum of the absolute terms divided by the absolute sum). Nothing
#' said above about the limits of the double range or about subnormal values
#' bounds the relative error of a `diff` whose terms cancel. For `diff`, a term
#' is zero when either
#' factor in `w * (A - B) * (A + B - 2 * h)` is exactly zero. As a fourth exception,
#' when the positive and negative terms of `diff`, scaled by the largest absolute
#' term, cancel to exactly 0, the function returns 0. A `diff` of 0 means the two
#' risks are not distinguished in double precision, not that they are equal.
#' A nonzero `diff` whose logarithm of the absolute value, computed in the same
#' way, is below `log(2^-1074)` stops the call with `representable` and the
#' field name.
