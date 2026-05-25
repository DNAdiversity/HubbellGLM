# Invert the mean of the Dirichlet process distribution

Given a target mean \\\mu\\ and a community size \\n\\, finds the
concentration parameter \\\alpha\\ such that \\E\[Y^{(n)}; \alpha\] =
\mu\\ via `uniroot`.

## Usage

``` r
inv_mean_dirichlet_process(mu_target, size)
```

## Arguments

- mu_target:

  Target mean value (scalar or vector, must be in \\(1, n)\\).

- size:

  Community size (positive integer, scalar or vector).

## Value

The concentration parameter \\\alpha\\ (scalar or vector).
