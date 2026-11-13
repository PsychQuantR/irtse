# irtse

**Status: unreleased.** `irtse` is a working name; the name, the interface and
the return values of `se_select()` may still change before a first release.

`irtse` compares two standard errors of an IRT ability estimate for dichotomous
logistic models with one or two item parameters (1P: all discriminations equal;
2P: discriminations may differ):

- the information-based standard error, `se_info()`;
- the parametric-bootstrap standard error, computed exactly by enumerating the
  law of the sufficient statistic rather than by Monte Carlo, `se_boot()`.

Two estimators are supported: the maximum likelihood estimate on responses that
are neither all correct nor all wrong (`"MLE"`), and Warm's weighted likelihood
estimate on all responses (`"WLE"`).

Exact standard errors are not new: `catR::semTheta(sem.exact = TRUE)` already
computes them. What `irtse` adds is

- `design_cumulants()`, `delta2()`, `delta3()`: closed forms of the second- and
  third-order coefficients of the difference between the scaled mean squared
  errors of the two standard errors;
- `risk_diff()`: the exact scaled mean squared errors at a given test length and
  ability, and their difference;
- `se_select()`: a provisional rule that chooses between the two standard
  errors from the signs of these coefficients at the ability estimate. The rule
  carries no guarantee at a fixed test length or near a zero of the
  second-order coefficient, and its `margin` argument has no theoretical
  default.

## Installation

Not yet published. `irtse` depends on the package `irtcore`.

## Exact enumeration

When all discriminations are equal, or are integer multiples of a common unit,
the law of the sufficient statistic is computed on a lattice. Otherwise every
response pattern is enumerated. Designs beyond the documented limits stop with
an error; the package never falls back to an approximation.
