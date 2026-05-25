# Variance-covariance matrix for a HubbellGLM fit, optionally adjusted for shared observations via a sandwich estimator.

Variance-covariance matrix for a HubbellGLM fit, optionally adjusted for
shared observations via a sandwich estimator.

## Usage

``` r
vcov_shared(fit, similarity = NULL)
```

## Arguments

- fit:

  An object of class `HubbellGLM`.

- similarity:

  An optional \\n \times n\\ similarity matrix encoding dependence
  between the \\n\\ observations used to fit `fit`. The diagonal must be
  1 and all off-diagonal entries must be in \\\[0, 1)\\. If `NULL`
  (default), returns `vcov(fit)`.

## Value

A \\p \times p\\ variance-covariance matrix.
