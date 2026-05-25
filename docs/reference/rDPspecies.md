# Simulate the number of distinct species from a Dirichlet process

Draws the number of distinct species \\Y^{(n)}\\ from the Dirichlet
process prior with concentration parameter \\\alpha\\ and sample size
\\n\\.

## Usage

``` r
rDPspecies(alpha, size)
```

## Arguments

- alpha:

  Concentration parameter (positive scalar or vector).

- size:

  Sample size \\n\\ (positive integer, scalar or vector).

## Value

An integer vector of simulated species counts.
