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
- `se_diagnostic()`: the same coefficients at the ability estimate, their sign,
  and the nearest zero of the second-order coefficient, without choosing a
  standard error;
- `se_select()`: a provisional rule that chooses between the two standard
  errors from the signs of these coefficients at the ability estimate. It is
  not a recommendation: the authors have not decided whether a switching rule
  should be used at all.

## What is known about the selection rule

For logistic designs and the MLE on truncated responses (the WLE part of the
underlying mathematical note has not been refereed):

1. On a set of abilities that stays away from the zeros of the second-order
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
4. Nothing is guaranteed at a fixed test length, and every statement assumes
   design averages that do not change with J (identical items, or a fixed block
   of items repeated).

## Installation

Not yet published. `irtse` depends on the package `irtcore`.

## Exact enumeration

When all discriminations are equal, or are integer multiples of a common unit,
the law of the sufficient statistic is computed on a lattice. Otherwise every
response pattern is enumerated. Designs beyond the documented limits stop with
an error; the package never falls back to an approximation.
