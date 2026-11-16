# irtse

**Status: unreleased.** The interface may still change before a first release.
The package is not on CRAN or any other repository.

`irtse` compares two standard errors of an IRT ability estimate for dichotomous
logistic models, `P_j(theta) = 1 / (1 + exp(-a_j (theta - b_j)))`, with one or
two item parameters (1P: all discriminations equal; 2P: discriminations may
differ). Two estimators are supported:

- `"MLE"`: the maximum likelihood estimate on responses that are neither all
  correct nor all wrong; every law used for the MLE is conditional on such
  responses;
- `"WLE"`: Warm's weighted likelihood estimate on all responses.

Exact standard errors are not new: `catR::semTheta(sem.exact = TRUE)` already
computes them. What `irtse` adds is the coefficients of the difference between
the two standard errors' scaled mean squared errors, the exact risks, and a
diagnostic built on them.

## What is exact and what is asymptotic

- **Exact computations** (no Monte Carlo, no approximation): `se_boot()` and
  `risk_diff()`. They enumerate the law of the sufficient statistic
  `T = sum(a * X)`: on a lattice when the discriminations are equal or are
  integer multiples of a common unit, and pattern by pattern otherwise. Item
  probabilities are computed on the log scale with `plogis(..., log.p = TRUE)`, so
  the law stays accurate at extreme abilities. `se_boot()` stops with an
  error above 100000 atoms, `risk_diff()` above 2000 atoms, and either above a
  lattice width of 1e6; the error states which limit and why. The package never
  falls back to an approximation.
- **Closed forms**: `se_info()` (the information-based standard error) and
  `design_cumulants()`.
- **Asymptotic coefficients**: `delta2()` and `delta3()` are the coefficients
  Delta_2 and Delta_3 of `J^(-2)` and `J^(-3)` in the expansion of the
  difference of the two scaled mean squared errors (the `J^(-1)` coefficient is
  zero). They describe long tests whose design averages do not change with J
  (identical items, or a fixed block of items repeated). The MLE and the WLE
  have different Delta_2; no Delta_3 is available for the WLE (`delta3()`
  returns `NA` with a message).

## Numerical accuracy and limits

Every exported function requires `1e-3 <= a_j <= 1e3`, `abs(b_j) <= 1e3`
and `abs(theta) <= 1e3`, inclusive. The bounds are the doubles represented by
the R literals `1e-3` and `1e3`. This range binds the supplied arguments only.
Estimates and internal evaluation points may lie outside it and are neither
rejected nor clamped, including root intervals and the window, grid points and
`nearest_zero` of `se_diagnostic()`. For example, one item with `a = 1e-3` at
its difficulty has WLE estimates at `b - 1000 * log(3)` and
`b + 1000 * log(3)`. `risk_diff(..., margin = )` evaluates the selection rule
at these estimates; supplying an estimate outside the range as `theta` to
`se_select()` is rejected. Numeric arguments (`a`, `b`, `theta`, `margin`) and
`estimator` that carry a class are rejected before any class method is called;
plain numeric vectors with names are accepted and their names are ignored.

All computations are in double precision. Inputs within the supported range
remain subject to the following numerical limits.

- At the supplied ability, information and coefficients use the
  logits `a * (theta - b)` computed from the supplied inputs in double precision.
  The law of T uses these logits with the identified discriminations described
  below.
  `se_boot()`, `risk_diff()` and `se_select()` keep estimates as offsets from the
  midpoint of the difficulties, so a common shift of every `b` and of `theta`
  changes the results only through rounding of the shifted inputs.
  Diagnostic zero searches use offsets from the supplied `theta`.
- A design average (`k3`, ..., `k7`) is a sum over items and can lose accuracy
  to cancellation between items: its relative error is of the order of 1e-15
  times the sum of the absolute terms divided by the absolute sum. Items that
  cancel exactly (same discrimination, opposite logits) are removed first, so a
  symmetric design gives exactly 0. When the remaining terms sum to exactly 0,
  the function stops with an error containing "not resolved", except in the zero
  search of `se_diagnostic()`, where an unresolved average counts as a zero.
  Accuracy also depends on each item's cumulant polynomial near its zeros and
  on the subtraction `theta - b` of nearly equal numbers.
  `delta2()` and `delta3()` inherit these errors
  and lose further accuracy near their own zeros.
- Log scales and normalised cumulant ratios reduce intermediate overflow and
  underflow. A call can still stop with "representable" when an intermediate
  value, such as a scaled cumulant ratio, is not representable, even if the
  returned quantity would be. The four scaled risks (`smse_info`, `smse_boot`,
  `diff`, `smse_select`) are summed on the log scale: each term is formed as a
  whole from the logarithm of its weight and of its two differences, so the sum
  never forms a weight, a squared difference or a partial product on its own,
  and none of them can stop a call or remove a term by leaving the double range.
  Apart from requiring that the log weights and differences it receives are
  finite, the sum checks the range only on the logarithm of the whole sum and on
  the value rebuilt from it.
  Weighted score terms are formed on the log scale, including the
  discrimination. The function and its derivative use the same scale at each
  point for signs, exact-zero tests and Newton steps, so the location of a root
  does not depend on the double range, even when every score term is smaller
  than the smallest positive double.
- A value computed as outside the range of double precision stops with an error
  containing
  "representable"; the package does not return `Inf`, `NaN`, or 0 in place of
  a non-zero value, except for a probability in `p_select` below the smallest
  positive double, the lattice identification described below, a risk term
  whose report agrees with its target to within rounding, or exactly cancelling
  terms of `diff`. In the third case,
  the squared difference is zero; for example,
  `risk_diff(c(1, 1 + 2^-36), c(0, 0), 0, "MLE")$smse_boot` can be zero even
  though the mathematical risk is positive. A risk whose logarithm, as this
  summation computes it, is below `log(2^-1074)` stops the call with
  `representable` and the field name, even when ordinary arithmetic would round
  the value up to `2^-1074`; a risk whose value rebuilt from that logarithm is
  not finite, or is 0, stops the call in the same way. The value of a scaled
  risk is
  rebuilt from its logarithm. Rounding that logarithm alone contributes a relative
  error of the order of 1e-13 where the logarithm is near 700 in absolute value;
  errors already in the reports and log weights add to it, and so does rounding
  the result to a double. A returned risk in the subnormal range is a multiple
  of `2^-1074`, so rounding the result to a subnormal double can give a much
  larger
  relative error: in `risk_diff(rep(1000, 3), c(-1.44, 0, 1.44), -0.72 + 5e-9,
  "MLE")` the exact `smse_boot` is 1.6 times `2^-1074` and the call returned 2
  times `2^-1074`. A
  subnormal risk
  is not claimed to be within `2^-1074` of the exact value, and no general
  bound on the relative error of a returned risk is claimed.
  Both decisions follow the logarithm computed from the reports, the target and
  the log weights in double precision, not the exact risk, and the two can lie
  on different sides of `2^-1074` or of the largest double. At `2^-1074` this is
  not confined to risks close to the limit: when a report and the target nearly
  agree,
  their difference keeps few significant digits and the computed risk can be
  wrong by a large factor. In `risk_diff(rep(1000, 3), c(-2 * r, 0, 2 * r), -r +
  2^-52, "MLE")` with `r` set to `0.685794632511088` the report and the target
  agree to
  a relative 1.1e-13, the exact `smse_boot` is 2.05 times `2^-1074`, which is
  representable, and the call stopped. With `r` set to `0.6907256300571373` and
  the ability `-r + 2^-49` the exact
  value is 0.96 times `2^-1074`, which is below the smallest positive double,
  and the call returned `2^-1074`. At the largest double, in `risk_diff(rep(512,
  2), c(-B, B), 0, "MLE")` with `B` set to `1042 * log(2) / 512` the exact
  `smse_info` is above the largest double by a relative 9.9e-15 and the call
  returned 1.7976931348622732e308. In general the package does not claim that
  the outcome agrees with the side of either limit on which the exact risk lies.
  The outcomes of the calls named here were measured with R 4.6.1 on aarch64
  macOS and can differ on another platform.
  For `diff`, cancellation multiplies the errors of its terms by its condition
  number (the sum of the absolute terms divided by the absolute sum). Nothing
  said above about the limits of the double range or about subnormal values
  bounds the relative error of a `diff` whose terms cancel. For `diff`, a term
  is zero when either factor
  in `w * (A - B) * (A + B - 2 * h)` is exactly zero. As a fourth exception, when
  the positive and negative terms of `diff`, scaled by the largest absolute term,
  cancel to exactly 0, the function returns 0. A `diff` of 0 means the two risks
  are not distinguished in double precision, not that they are equal.
  A nonzero `diff` whose logarithm of the absolute value, computed in the same
  way, is below `log(2^-1074)` stops the call with `representable` and the field
  name.
  For example, `delta2(1, 0, 300, "MLE")` stops because the coefficient is not
  representable. `delta2(1, 0, 2^-1074, "MLE")` also stops with "representable"
  because `tanh(eta / 2)` is zero for a nonzero logit. A non-finite logit also
  causes this error, though it cannot arise at the supplied ability within
  the supported range, where logits have absolute value at most 2e6.
- Discriminations are identified with a commensurable design when multiplying
  their ratios to the smallest discrimination by an integer from 1 to 1000
  puts them within eight relative machine epsilons of integers. This includes
  decimal inputs such as `c(0.1, 0.2, 0.3)` and distinct supplied doubles such
  as `c(1, 1 + 2^-50)`; otherwise every response pattern is its own atom.
  The identified design is used throughout the law-based functions: totals,
  estimating equations, estimates, the law at `theta` and at each estimate,
  and every report and rule coefficient at an estimate. The `se_info`, `delta2`
  and `delta3` fields of `se_select()` still use the supplied discriminations.
  Identified totals form one atom and can give an exactly zero bootstrap SE.
- The WLE equation can have several roots. The package counts roots on a grid
  together with the extrema of the function that the grid reveals, and stops
  with an error naming T when the count is not one. The count is certified
  only to the resolution of the grid.

## Usage

```r
library(irtse)

# Six identical Rasch items at their difficulty
a <- rep(1, 6)
b <- rep(0, 6)

# The two standard errors at an ability estimate (actual SEs, not scaled)
se_info(a, b, theta = 0)                     # 0.8165
se_boot(a, b, theta = 0, estimator = "MLE")  # 0.8566

# Design averages and the risk coefficients
design_cumulants(a, b, theta = 0)            # k2 ... k7
delta2(a, b, theta = qlogis(3/10), estimator = "MLE")  # -48800/3087
delta3(a, b, theta = 0, estimator = "MLE")             # -121/3

# Exact scaled risks at a test length and true ability
r <- risk_diff(a, b, theta = 0, estimator = "MLE")
r$smse_info  # 0.0696
r$smse_boot  # 0.0972
r$diff       # smse_info - smse_boot, not multiplied by J^2

# Diagnostic without choosing: coefficients, their sign, nearest zero of Delta_2
se_diagnostic(rep(c(1, 2), 25), rep(c(-1/2, 1), 25), theta = 0.03,
              estimator = "MLE")
```

`se_info()`, `se_boot()` and `se_select()` return actual standard errors;
`risk_diff()` works with the scaled reports, multiplied by `sqrt(J)`.

## The selection rule is provisional

> `se_select()` chooses between the two standard errors from the sign of
> Delta_2 at the estimate (and of Delta_3 for the MLE within `margin`). It is
> not a recommendation, its form is provisional, and `margin` has no
> theoretical default. The author has not yet decided how the selection rule
> will be presented, and whether a switching rule will be used at all.
> `risk_diff(..., margin = )` returns the exact scaled risk of the rule
> itself, so its behaviour can be checked for a given design.

What is known, for logistic designs and the MLE on truncated responses (the
WLE part of the underlying mathematical note has not been refereed):

1. On a compact set of abilities that contains no zero of the second-order
   coefficient Delta_2, choosing by the sign of Delta_2 at the estimate gives a
   scaled risk within an exponentially small amount of the smaller of the two
   fixed scaled risks.
2. Any rule that switches between the two standard errors on the estimate has,
   near the abilities where it switches, a scaled risk equal to the smaller
   fixed risk plus or minus a term of order J^(-3/2). The sign of that term is
   the sign of c G' Delta_2' at the switch point and cannot be chosen by the
   user. Always reporting one standard error costs at most |D_J|, of order
   J^(-2). A margin moves the switch to the edge of the margin band; it does not
   remove the term.
3. Example: a 2P design with half the items (a, b) = (1, -1/2) and half
   (2, 1), at ability 0.01805334808, near a zero of Delta_2. Choosing by the
   sign of Delta_2 is worse than both fixed standard errors, by 18.1%, 16.9%,
   9.1% and 4.6% of the smaller fixed scaled risk at J = 50, 100, 400 and 1600.
   At J = 50, `risk_diff(..., margin = 0)` gives 0.0382 (information),
   0.0386 (bootstrap) and 0.0451 (the rule).
4. Nothing is guaranteed at a fixed test length, and every statement assumes
   design averages that do not change with J (identical items, or a fixed block
   of items repeated).

## Installation

Not yet published. `irtse` depends only on base R (`stats`).

## Tests

From the source folder: `Rscript -e 'devtools::test(".")'`. The expected
values in the tests are exact rationals and certified values from the
mathematical notes behind the package.
