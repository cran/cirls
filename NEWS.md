# 0.3.0

## New features
- Added `check_cmat` and `coef_simu` to the list of exported functions as they can be useful for specific use cases.
- Added full documentation.

## Bug fixes
- In `check_cmat`, removed the call to `limSolve::nnls()` to be replaced by `coneproj::coneB` (also a NNLS solver) to reduce the number of dependencies.

# 0.2.1

## New features
- Initialization of a short documentation for several functions.

## Bug fixes
- `cirls.control` now checks for constraint matrix irreducibility.
- `cirls.control` is now exported.
- changed the function to determine redundant constraints in `check_cmat`
- No warning message on row rank of `Cmat` anymore
- `Cmat` is now checked only when there are more than one row

# 0.2.0

## New features
- Method vcov for cirls to compute corrected covariance matrices
- Method confint for cirls to compute feasible confidence intervals
- Added several QP solvers: quadprog (the original one), osqp and coneproj.

## Changes
- A warning is now displayed when Cmat is not of full row rank
- vcov and confint return NA matrices if Cmat is not of full row rank
- Changed residual df computation to account for active constraints
- Replaced bvec by lb (lower bound) and ub (upper bound). Allows equality constraints.
- Added cirls class to glm output

## Bug fixes
- cirls.fit has the same behaviour as glm.fit when model is singular: fill with NAs
