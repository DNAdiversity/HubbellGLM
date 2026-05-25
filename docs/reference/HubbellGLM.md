# Hubbell regression This function mimics the GLM structure. We do so to avoid conflicts with the current glm implementation (if it gets updated)

Hubbell regression This function mimics the GLM structure. We do so to
avoid conflicts with the current glm implementation (if it gets updated)

## Usage

``` r
HubbellGLM(
  formula,
  family = hubbell(sigma = 0),
  data,
  weights,
  subset,
  na.action,
  start = NULL,
  etastart,
  mustart,
  offset,
  control = list(...),
  model = TRUE,
  x = FALSE,
  y = TRUE,
  singular.ok = TRUE,
  contrasts = NULL,
  ...
)
```

## Arguments

- formula:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- family:

  a member of the hubbell family for `glm.fit.hubbell`.

- data:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- weights:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- subset:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- na.action:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- start:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- etastart:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- mustart:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- offset:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- control:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- model:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- x:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- y:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- singular.ok:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- contrasts:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

- ...:

  Same as [`stats::glm`](https://rdrr.io/r/stats/glm.html)

## Value

An object of class `c('HubbellGLM', 'lm')`, which inherits from the
class "lm"
