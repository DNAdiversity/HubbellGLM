# Barro Colorado Island subplot data

A dataset of 50 subplots from the Barro Colorado Island (BCI) 50-ha
permanent forest plot in Panama, one of the most intensively studied
tropical forests in the world.

## Usage

``` r
BarroColorado
```

## Format

A `data.frame` with 50 rows and 5 variables:

- n:

  Community size: total number of individual trees in the subplot
  (positive integer).

- y:

  Observed species richness: number of distinct species in the subplot
  (positive integer, \\1 \le y \le n\\).

- EnvHet:

  Environmental heterogeneity index (continuous, higher values indicate
  greater habitat complexity).

- Habitat:

  Habitat type (factor with 5 levels: `OldHigh`, `OldLow`, `OldSlope`,
  `Swamp`, `Young`).

- Stream:

  Proximity to a stream (`"Yes"` / `"No"`).

## Source

Condit et al. (2002). Beta-diversity in tropical forest trees.
*Science*, 295(5555), 666–669.
