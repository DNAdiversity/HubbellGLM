# Stepwise model selection for HubbellGLM models

Stepwise model selection for HubbellGLM models

## Usage

``` r
# S3 method for class 'HubbellGLM'
step(
  object,
  scope,
  scale = 0,
  direction = c("both", "backward", "forward"),
  trace = 1,
  keep = NULL,
  steps = 1000,
  k = 2,
  ...
)
```

## Arguments

- object:

  A fitted `HubbellGLM` object.

- scope:

  A formula or list of formulas defining the range of models to
  consider. See [`step`](https://rdrr.io/r/stats/step.html) for details.

- scale:

  Scaling factor for the AIC penalty. Default is 0.

- direction:

  Direction of stepwise search: `"both"`, `"backward"`, or `"forward"`.

- trace:

  Integer; if positive, prints information at each step.

- keep:

  A function to retain extra information at each step.

- steps:

  Maximum number of steps. Default is 1000.

- k:

  Penalty multiplier for the number of parameters (2 for AIC).

- ...:

  Additional arguments passed to the fitting function.
