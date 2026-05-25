# Function to predict the response from a HubbellGLM model.

Function to predict the response from a HubbellGLM model.

## Usage

``` r
# S3 method for class 'HubbellGLM'
predict(
  object,
  newdata = NULL,
  type = c("link", "response", "terms"),
  se.fit = FALSE,
  dispersion = NULL,
  terms = NULL,
  na.action = na.pass,
  ...
)
```
