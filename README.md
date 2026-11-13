# irtse

**Status: unreleased.** `irtse` is a working name; the name, the interface and
the return values of `se_select()` may still change before a first release.
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
  probabilities are computed as `P = plogis(eta)` and `Q = plogis(-eta)`, so
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
