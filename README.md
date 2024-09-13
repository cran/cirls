# cirls: Constrained Iteratively Reweighted Least Squares

<!-- badges: start -->
  [![R-CMD-check](https://github.com/PierreMasselot/cirls/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PierreMasselot/cirls/actions/workflows/R-CMD-check.yaml)
 <!-- badges: end -->

The package `cirls` provides routines to fit Generalized Linear Models (GLM) with coefficients subject to linear constraints, through a constrained iteratively reweighted least-squares algorithm. 

## Installation

The easiest way to install the `cirls` package is to install it from CRAN

```R

install.packages("cirls")

```

Although not recommended, the development version can be installed from GitHub using the `devtools` package as

```R
devtools::install_github("PierreMasselot/cirls")
```

## Usage

The central function of the package is `cirls.fit` meant to be passed through the `method` argument of the `glm` function. The user is also expected to pass a either constraint matrix or a list of constraint matrices through the `Cmat` argument, and optionally lower and upper bound vectors `lb` and `ub`. 

The package also contains dedicated methods to extract the variance-covariance matrix of the coefficients `vcov. cirls` as well as confidence intervals `confint.cirls`.

The example below show how to use the package to perform nonnegative regression. See `?cirls.fit` for more comprehensive examples.

```R
# Simulate predictors and response with some negative coefficients
set.seed(111)
n <- 100
p <- 10
betas <- rep_len(c(1, -1), p)
x <- matrix(rnorm(n * p), nrow = n)
y <- x %*% betas + rnorm(n)

# Define constraint matrix
Cmat <- diag(p)

# Fit GLM by CIRLS
res <- glm(y ~ x, method = cirls.fit, Cmat = list(x = Cmat))
coef(res)

# Obtain vcov and confidence intervals
vcov(res)
confint(res)
```

## References

*To come*
