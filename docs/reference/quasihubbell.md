# Quasi-Hubbell family of generalized linear models

Quasi-Hubbell family of generalized linear models

## Usage

``` r
quasihubbell(link = "polyseries", sigma = 0)
```

## Arguments

- link:

  The link function to use. Available options are `dp`, `py` and
  `polyseries`

- sigma:

  Hyperparameter determining the polynomial growth. Must be `sigma < 1`
  for link `link = 'polyseries'`. It is set to 0 when `link = 'dp'`
  (logarithmic growth).
