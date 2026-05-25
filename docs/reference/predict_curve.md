# Predict the accumulation curve from a HubbellGLM model for new observations.

Predict the accumulation curve from a HubbellGLM model for new
observations.

## Usage

``` r
predict_curve(fit, xnew, n = 3000, npoints = 100, .vcov = NULL)
```

## Arguments

- fit:

  An object of class `HubbellGLM`.

- xnew:

  A `data.frame` with the same predictor columns used to fit `fit`. Each
  row produces one set of accumulation curve predictions.

- n:

  Total number of individuals to extrapolate to. Must be a single
  positive integer.

- npoints:

  Number of equally-spaced grid points between 1 and `n`. Must satisfy
  `npoints <= n`. Default is 100.

- .vcov:

  Optional variance-covariance matrix for the regression coefficients.
  Must be a square matrix of dimension equal to the number of
  coefficients in `fit`. Defaults to `vcov(fit)`.

## Value

A named list with three elements:

- `n`:

  Integer vector of length `npoints` with the grid points (number of
  individuals).

- `mean`:

  A `data.frame` with `npoints` rows and one column per row of `xnew`,
  containing the predicted species richness at each grid point.

- `se`:

  A `data.frame` of the same shape as `mean`, containing standard errors
  computed via the delta method.
