# Summary function for the HubbellGLM method

Summary function for the HubbellGLM method

## Usage

``` r
# S3 method for class 'HubbellGLM'
summary(
  object,
  dispersion = NULL,
  correlation = FALSE,
  symbolic.cor = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted `HubbellGLM` object.

- dispersion:

  Optional dispersion parameter. If `NULL`, estimated from the data.

- correlation:

  Logical; if `TRUE`, print the correlation matrix of the estimated
  coefficients. Default is `FALSE`.

- symbolic.cor:

  Logical; if `TRUE`, print correlations as symbols. Default is `FALSE`.

- ...:

  Additional arguments (currently unused).
