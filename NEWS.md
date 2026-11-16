# irtse 0.0.0.9000

Unreleased development version.

* Correct what the documents say about the limits of the double range. A scaled
  risk stops the call when the logarithm the summation computes is below
  `log(2^-1074)` or when the value rebuilt from it is not finite or is 0. That
  logarithm
  is computed from rounded reports, target and log weights, so the outcome need
  not agree with the side of the limit on which the exact risk lies; when a
  report and the target nearly agree the computed risk can be wrong by a large
  factor. The earlier text claimed that this happens only within a relative
  1e-12 of either limit; no width is claimed now. No function changed.

* Clarify scaled-risk accuracy near both ends of the double range: rounding the
  logarithm from which a risk is rebuilt contributes a relative error of the
  order of 1e-13, rounding the result to a subnormal double can contribute much
  more, and subnormal risks are not guaranteed to be within one spacing.
  Cancellation in `diff` amplifies the errors of its terms.

* Sum all four scaled risks on the log scale, retaining representable terms
  despite underflowing weights or overflowing squared reports. Stop when the
  logarithm computed for a risk is below `log(2^-1074)`, even when ordinary
  arithmetic would round the value up, and name the field in range errors.
  Exactly cancelling `diff` terms return 0: the two
  risks are not distinguished in double precision.

* Reject classed numeric arguments and classed `estimator` values before any
  class method is called; continue to accept named plain numeric vectors.

* Form weighted score terms on the log scale and locate roots from a common
  scale for the function and derivative, including beyond the double range.

* Stop for scaled risks whose computed value underflows, including the selection
  risk, while preserving exact zeros from rounded reports or cancelling
  differences.

* Reject numeric arguments carrying a class before arithmetic; continue to
  accept named plain numeric vectors.

* Retain score tails representable only as subnormal numbers by obtaining
  probabilities from their logarithms, correcting MLE roots and bootstrap SEs
  for separated difficulties.

* Use the identified discriminations throughout law-based computations,
  including estimating equations, laws, reports and rule coefficients at
  estimates. Keep supplied discriminations for information and coefficients
  at the supplied ability.

* Restrict supplied inputs to 1e-3 <= a_j <= 1e3, |b_j| <= 1e3 and |theta| <= 1e3 (inclusive); estimates and internal evaluation points remain unrestricted by these bounds.

* Preserve supplied-ability logits, search diagnostic zeros in offsets from theta, reject unrepresentable logits, and clarify accuracy and lattice identification.

* Numerical range and root accuracy: retain small-logit odd cumulants and
  saturated score tails, compute information independently of higher cumulants,
  and use log scales for extreme abilities and discriminations. Law-based
  calculations keep estimates in shifted design coordinates.
  Exact opposite item pairs cancel before summation; unresolved cancellation
  stops with "not resolved". WLE root counts include function extrema found
  through derivative sign changes and searches at sampled derivative extrema.
  Inputs within eight relative machine epsilons of a commensurable design are
  identified with it, including distinct supplied doubles. Selection
  and its exact risk use the same coefficients, and a no-choice atom makes the
  selection risk undefined even when its probability underflows. Documentation
  states the numerical and root-search limits and the compact-set scope of the
  provisional rule.

* `se_info()`: information-based standard error at an ability estimate.
* `se_boot()`: exact parametric-bootstrap standard error, by exact enumeration
  of the law of the sufficient statistic (lattice convolution for equal or
  commensurable discriminations, pattern enumeration otherwise; `P` and `Q`
  both computed with `plogis(..., log.p = TRUE)`); MLE conditional on responses that are neither
  all correct nor all wrong, WLE on all responses.
* `design_cumulants()`: design averages k2 to k7.
* `delta2()`, `delta3()`: closed forms of the second- and third-order
  coefficients of the difference of the two scaled mean squared errors; the
  MLE and the WLE have separate Delta_2, and `delta3()` returns `NA` for the
  WLE.
* `risk_diff()`: exact finite-length scaled mean squared errors and their
  difference; with `margin`, also the exact scaled risk of the selection rule.
* `se_diagnostic()`: Delta_2, Delta_3, the sign of Delta_2 and the nearest zero
  of Delta_2, without choosing a standard error.
* `se_select()`: provisional selection rule with a required `margin`; not a
  recommendation (see the README for what is and is not known about it).
* Exact enumeration stops with an error above its limits (100000 atoms for
  `se_boot()`, 2000 for `risk_diff()`, lattice width 1e6), stating the limit
  and the reason, instead of approximating.
* No dependency other than `stats`.
