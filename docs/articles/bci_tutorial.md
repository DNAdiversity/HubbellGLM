# Getting started with HubbellGLM: the Barro Colorado Island dataset

### Introduction

The **HubbellGLM** package implements Hubbell regression, a generalized
linear model (GLM) for $\alpha$-diversity (Zito et al. 2026). The model
links environmental covariates to Hubbell’s fundamental biodiversity
number $\alpha$, a single dimension-free parameter that controls the
rate at which new species accumulate as a function of sample size, and
can be converted analytically into Shannon entropy, the Simpson index,
and Hill numbers. This bypassess the need to choose which biodiversity
index to use in the analysis.

This vignette introduces the statistical framework and demonstrates the
package using the **Barro Colorado Island (BCI) dataset**, a classic
tree-census dataset from a 50-ha tropical forest plot in Panama. The
original version of the data is downloaded from the R package `vegan`
(Oksanen et al. 2024).

------------------------------------------------------------------------

### Model details

#### Hubbell’s fundamental biodiversity number

Let $Y^{(n)}$ denote the number of distinct species (richness) observed
in a sample of $n$ individuals drawn from a community. Under Hubbell’s
neutral theory of biodiversity (Hubbell 2001), the distribution of
$Y^{(n)}$ is given by

$${\mathbb{P}}\left( Y^{(n)} = y;\,\alpha \right) = \frac{\alpha^{y}\,\Gamma(\alpha)}{\Gamma(\alpha + n)}\,\left| s(n,y) \right|,\qquad y = 1,\ldots,n,$$

where $\left| s(n,y) \right|$ are the (unsigned) Stirling numbers of the
first kind, $\Gamma(\cdot)$ is the gamma function, and $\alpha > 0$ is
the **fundamental biodiversity number** (Hubbell 2001). This
distribution arises from the Ewens sampling formula (Ewens 1972), which
is mathematically equivalent to the Dirichlet process described in
(Ferguson 1973) and (Antoniak 1974).

In this model, the richness $Y^{(n)}$ accumulates logarithmically as a
function of $n$, so that
$$\lim\limits_{n\rightarrow\infty}\frac{Y^{(n)}}{\log n} = \alpha.$$
This relationship makes $\alpha$ equivalent to the biodiveristy number
in the log-series model proposed by (Fisher, Corbet, and Williams 1943).
A larger $\alpha$ means richer, more even communities; a smaller
$\alpha$ means communities dominated by a few abundant species. In
particular, the expected species richness in a sample of size $n$ is

$$\mu = {\mathbb{E}}\left( Y^{(n)} \right) = \sum\limits_{j = 0}^{n - 1}\frac{\alpha}{\alpha + j}.$$

Hence, the above mean is interpreted as a **species accumulation curve**
that grows logarithmically with $n$. Moreover, since Hubbell’s model is
equivalent to the Dirichlet process, all standard biodiversity indices
have closed-form model-based expressions in terms of $\alpha$(Rigon,
Hsu, and Dunson 2025; Pitman 2006):

| Index | Sample-based | Model-based |
|:---|:--:|:--:|
| $\alpha$-diversity | ${\widehat{\alpha}}_{F}\log\left( 1 + n/{\widehat{\alpha}}_{F} \right) = y$ | $y = \sum_{j = 0}^{n - 1}\frac{\alpha}{\alpha + j}$ |
| Shannon | $-\sum_{k}p_{k}\log p_{k}$ | $\psi(\alpha + 1) - \psi(1)$ |
| Simpson | $\sum_{k}p_{k}^{2}$ | $1/(\alpha + 1)$ |
| Hill ($q$) | $\frac{1}{q - 1}\left( 1 - \sum_{k}p_{k}^{q} \right)$ | $\alpha B(\alpha,q)$ |

Here $\psi(\cdot)$ is the digamma function, $B(\cdot,\cdot)$ is the beta
function, and $p_{k}$ is the relative abundance of the $k$-th species.

#### The Hubbell regression

Since the distribution of $Y^{(n)}$ is an exponential family, we can use
it to build a generalized linear model. In particular, our goal is to
link the fundamental biodiversity number $\alpha$ with
$\mathbf{x} \in {\mathbb{R}}^{p}$. We call this model **Hubbell
regression**. Specifically, letting $i = 1,\ldots,N$ denote the sampling
sites (i.e. the observations), we let

$$Y_{i}^{(n_{i})} \sim \text{Hubbell}\left( \alpha_{i},n_{i} \right),\qquad\log\alpha_{i} = \eta_{i} = \mathbf{x}_{i}^{\top}{\mathbf{β}}.$$

This is the **canonical link** ($\sigma = 0$). A one-unit increase in
covariate $x_{ip}$ translates into a
$100\%\left( e^{\beta_{p}} - 1 \right)$ change in $\alpha$-diversity.

##### The polynomial link function

The canonical link constrains species richness to grow logarithmically
with $n$. To allow for heavier-tailed abundance distributions
(e.g. arthropods), the package supports a more general **polynomial
link** controlled by $\sigma < 1$:

$$\mu_{i} = \sum\limits_{j = 0}^{n_{i} - 1}\frac{e^{\eta_{i}}}{e^{\eta_{i}} + j^{1 - \sigma}},\qquad\sigma < 1.$$

- $\sigma = 0$: logarithmic growth (canonical Hubbell link).
- $\sigma \in (0,1)$: polynomial growth; the **regression-based
  diversity index** is $S_{\sigma}(\eta) = e^{\eta}/\sigma$.
- $\sigma < 0$: finite asymptotic richness.

The parameter $\sigma$ is estimated by maximum likelihood via
[`estimate_sigma()`](https://alessandrozito.github.io/HubbellGLM/reference/estimate_sigma.md)
and then held fixed across nested model specifications to ensure
comparable regression coefficients.

------------------------------------------------------------------------

## Robust standard errors: Jaccard-adjusted variance–covariance

Standard GLM inference assumes conditionally independent observations.
In species-richness data, the same species may appear in multiple sites,
inducing cross-site dependence that the likelihood ignores. HubbellGLM
corrects for this using **heteroskedastic-consistent sandwich standard
errors**, where spatial dependence between sites $i$ and $k$ is measured
by the Jaccard similarity index $s_{ik}$ (the fraction of shared species
between the two samples).

Let $\mathbf{X} \in {\mathbb{R}}^{N \times p}$ be the design matrix,
$\widehat{\mathbf{W}} = \text{diag}\left( {\widehat{w}}_{1},\ldots,{\widehat{w}}_{N} \right)$
the matrix of IRLS working weights, and
${\widehat{\mathbf{u}}}_{i} = \partial\log\mathcal{L}({\mathbf{β}})/\partial{\mathbf{β}}|_{\widehat{\mathbf{β}}}$
the score for observation $i$. The **Jaccard-adjusted
variance–covariance matrix** is

$$\widehat{\text{se}}\left( \widehat{\mathbf{β}} \right) = \left( \mathbf{X}^{\top}\widehat{\mathbf{W}}\mathbf{X} \right)^{-1}\left( \sum\limits_{i,k}s_{ik}\,{\widehat{\mathbf{u}}}_{i}{\widehat{\mathbf{u}}}_{k}^{\top} \right)\left( \mathbf{X}^{\top}\widehat{\mathbf{W}}\mathbf{X} \right)^{-1}.$$

This is computed by `vcov_species(fit, Dshared)`, where `Dshared` is an
$N \times N$ matrix of pairwise Jaccard similarities. When species
membership between sites is unavailable (as in the BCI example below),
the standard Fisher-information-based `vcov(fit)` is used.

------------------------------------------------------------------------

## Package functions

### `HubbellGLM()` — fit a Hubbell regression

The main fitting function mirrors the syntax of
[`stats::glm`](https://rdrr.io/r/stats/glm.html). The response is
`cbind(n, y)`, where `n` is the sample size (number of individuals) and
`y` is the observed richness.

``` r
fit <- HubbellGLM(cbind(n, y) ~ x1 + x2,
                  family = hubbell(sigma = 0),
                  data   = mydata)
```

Key arguments:

| Argument | Description |
|:---|:---|
| `formula` | As in `glm`. Response must be `cbind(n, y)`. |
| `family` | A Hubbell family object from [`hubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/hubbell.md) or [`quasihubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/quasihubbell.md). |
| `data` | A data frame. |

The function returns an object of class `HubbellGLM` that is compatible
with [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`deviance()`](https://rdrr.io/r/stats/deviance.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html), and
[`step()`](https://rdrr.io/r/stats/step.html).

### `hubbell()` — family object

``` r
hubbell(link = "polyseries", sigma = 0)
```

Sets the GLM family. `sigma = 0` gives the canonical (logarithmic) link;
any value in $(-\infty,1)$ gives the polynomial link. Use
[`quasihubbell()`](https://alessandrozito.github.io/HubbellGLM/reference/quasihubbell.md)
for quasi-likelihood inference when the variance is suspected to be
inflated.

### `estimate_sigma()` — estimate the growth-rate parameter

``` r
sigma_hat <- estimate_sigma(cbind(n, y) ~ 1, data = mydata)
```

Finds the maximum-likelihood $\widehat{\sigma}$ by minimising the BIC of
the null model over $\sigma \in (-\infty,0.9)$ using `nlminb`. Fix
$\sigma$ to this value for all subsequent models to ensure comparability
of regression coefficients across nested specifications.

### `predict.HubbellGLM()` — predict richness at new sites

``` r
predict(fit, newdata = new_sites, type = "response")
```

- `type = "link"` returns
  $\widehat{\eta} = \mathbf{x}^{\top}\widehat{\mathbf{β}}$.
- `type = "response"` returns the predicted mean richness
  $\widehat{\mu} = g^{-1}\left( \widehat{\eta} \right)$.

### `predict_curve()` — predicted accumulation curve

``` r
curve <- predict_curve(fit, xnew = new_site, n = 2000, npoints = 200)
```

Returns a data frame with columns `n` (sample size grid), `pred`
(predicted richness), and `se` (standard error via the delta method).
Useful for visualising how richness accumulates at a specific new
location.

### `vcov_species()` — Jaccard-adjusted variance–covariance

``` r
Dshared  <- get_shared_species(species_matrix)  # N x N Jaccard matrix
vcov_adj <- vcov_species(fit, Dshared)
```

Returns the sandwich variance–covariance matrix described in Section 3.
Pass it to `lmtest::coeftest(fit, vcov = vcov_adj)` to obtain robust
z-scores and p-values.

------------------------------------------------------------------------

## BCI example

### The data

The `BarroColorado` dataset contains 50 subplots from the Barro Colorado
Island 50-ha plot. Each row is one subplot with:

- `n`: number of individual trees sampled.
- `y`: number of distinct species observed.
- `EnvHet`: environmental heterogeneity (continuous).
- `Habitat`: habitat type (factor with 5 levels).
- `Stream`: proximity to a stream (`Yes`/`No`).

``` r
library(HubbellGLM)

data("BarroColorado")
head(BarroColorado)
#>     n   y EnvHet  Habitat Stream
#> 1 448  93 0.6272 OldSlope    Yes
#> 2 435  84 0.3936   OldLow    Yes
#> 3 463  90 0.0000   OldLow     No
#> 4 508  94 0.0000   OldLow     No
#> 5 505 101 0.4608 OldSlope     No
#> 6 412  85 0.0768   OldLow     No
```

``` r
plot(BarroColorado$n, BarroColorado$y,
     col  = as.integer(BarroColorado$Habitat),
     pch  = 16, cex = 0.9,
     xlab = "Sample size (n)",
     ylab = "Observed richness (y)",
     main = "BCI: richness vs. sample size")
legend("bottomright", legend = levels(BarroColorado$Habitat),
       col = 1:5, pch = 16, cex = 0.8)
```

![Observed species accumulation: richness vs. sample size across habitat
types.](bci_tutorial_files/figure-html/eda-1.png)

Observed species accumulation: richness vs. sample size across habitat
types.

### The canonical link (logarithmic growth)

These are three models run with the canonical link

``` r
# M0: the null model
fit_null <- HubbellGLM(cbind(n, y) ~ 1,
                     family = hubbell(sigma = 0),
                     data   = BarroColorado)

# M1: environmental heterogeneity only
fit_M1 <- HubbellGLM(cbind(n, y) ~ EnvHet,
                     family = hubbell(sigma = 0),
                     data   = BarroColorado)

# M2: add habitat type
fit_M2 <- HubbellGLM(cbind(n, y) ~ EnvHet + Habitat,
                     family = hubbell(sigma = 0),
                     data   = BarroColorado)

# M3: add stream proximity
fit_M3 <- HubbellGLM(cbind(n, y) ~ EnvHet + Habitat + Stream,
                     family = hubbell(sigma = 0),
                     data   = BarroColorado)

summary(fit_M3)
#> 
#> Call:
#> HubbellGLM(formula = cbind(n, y) ~ EnvHet + Habitat + Stream, 
#>     family = hubbell(sigma = 0), data = BarroColorado)
#> 
#> Coefficients:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)      3.43326    0.05523  62.161  < 2e-16 ***
#> EnvHet           0.07146    0.09230   0.774  0.43880    
#> HabitatOldLow    0.14270    0.05534   2.579  0.00992 ** 
#> HabitatOldSlope  0.13352    0.06224   2.145  0.03193 *  
#> HabitatSwamp     0.26655    0.11023   2.418  0.01560 *  
#> HabitatYoung    -0.04071    0.10360  -0.393  0.69439    
#> StreamYes       -0.10782    0.05879  -1.834  0.06665 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for hubbell family taken to be 1)
#> 
#>     Null deviance: 49.331  on 49  degrees of freedom
#> Residual deviance: 34.543  on 43  degrees of freedom
#> AIC: 343.26
#> 
#> Number of Fisher Scoring iterations: 3
```

#### Interpreting coefficients

Coefficients are on the log-$\alpha$ scale. A coefficient
${\widehat{\beta}}_{p}$ means a one-unit increase in $x_{p}$ is
associated with a $100\%\left( e^{{\widehat{\beta}}_{p}} - 1 \right)$
change in $\alpha$-diversity.

``` r
round(100 * (exp(coef(fit_M3)) - 1), 2)
#>     (Intercept)          EnvHet   HabitatOldLow HabitatOldSlope    HabitatSwamp 
#>         2997.75            7.41           15.34           14.28           30.55 
#>    HabitatYoung       StreamYes 
#>           -3.99          -10.22
```

#### Model comparison via deviance

``` r
cat("Null deviance:    ", round(deviance(fit_null), 1), "\n")
#> Null deviance:     49.3
cat("M1 deviance:      ", round(deviance(fit_M1),  1), "\n")
#> M1 deviance:       49.2
cat("M2 deviance:      ", round(deviance(fit_M2),  1), "\n")
#> M2 deviance:       37.9
cat("M3 deviance:      ", round(deviance(fit_M3),  1), "\n")
#> M3 deviance:       34.5
```

The pseudo-$R^{2}$ for each model relative to the null:

``` r
1 - deviance(fit_M3) / deviance(fit_null)
#> [1] 0.2997722
```

#### Predicted accumulation curves

[`predict_curve()`](https://alessandrozito.github.io/HubbellGLM/reference/predict_curve.md)
returns the fitted accumulation curve (with pointwise standard errors)
at any new covariate combination. Here we compare a high-heterogeneity
`OldHigh` subplot with a near-stream `Young` subplot.

``` r
new_sites <- data.frame(
  n       = c(500, 500),
  y       = c(1L,  1L),      # placeholder required by the formula
  EnvHet  = c(0.6,  0.0),
  Habitat = factor(c("OldHigh", "OldSlope"),
                   levels = levels(BarroColorado$Habitat)),
  Stream  = factor(c("No", "Yes"),
                   levels = levels(BarroColorado$Stream))
)

curve1 <- predict_curve(fit_M3, xnew = new_sites[1, ], n = 600, npoints = 100)
curve2 <- predict_curve(fit_M3, xnew = new_sites[2, ], n = 600, npoints = 100)

pred1 <- curve1$mean[, 1]; se1 <- curve1$se[, 1]
pred2 <- curve2$mean[, 1]; se2 <- curve2$se[, 1]

plot(curve1$n, pred1, type = "l", col = "steelblue", lwd = 2,
     ylim = range(c(pred1 - 2*se1, pred2 + 2*se2)),
     xlab = "Sample size (n)", ylab = "Predicted richness",
     main = "Predicted accumulation curves")
lines(curve2$n, pred2, col = "tomato", lwd = 2)

# 95% CI bands
polygon(c(curve1$n, rev(curve1$n)),
        c(pred1 - 2*se1, rev(pred1 + 2*se1)),
        col = adjustcolor("steelblue", 0.15), border = NA)
polygon(c(curve2$n, rev(curve2$n)),
        c(pred2 - 2*se2, rev(pred2 + 2*se2)),
        col = adjustcolor("tomato", 0.15), border = NA)

legend("bottomright",
       legend = c("OldHigh, high EnvHet, no stream",
                  "Young, low EnvHet, near stream (OldSlope fixed)"),
       col = c("steelblue", "tomato"), lwd = 2, cex = 0.8)
```

![Predicted accumulation curves for two contrasting subplot
types.](bci_tutorial_files/figure-html/pred-curves-1.png)

Predicted accumulation curves for two contrasting subplot types.

#### Derived biodiversity indices

Once $\widehat{\alpha}$ is recovered from the fitted model, all standard
biodiversity indices follow analytically.

``` r
alpha_n <- function(mu, n) {
  uniroot(function(a) mean_dirichlet_process(a, n) - mu,
          c(1e-6, 1e6))$root
}

# Predictions at n = 500 for each subplot
mu_pred <- predict(fit_M3, type = "response")
n_vec   <- BarroColorado$n

alpha_vals   <- mapply(alpha_n, mu_pred, n_vec)
shannon_vals <- digamma(alpha_vals + 1) - digamma(1)
simpson_vals <- 1 / (alpha_vals + 1)

cat("Median alpha:   ", round(median(alpha_vals),   3), "\n")
#> Median alpha:    35.729
cat("Median Shannon: ", round(median(shannon_vals), 3), "\n")
#> Median Shannon:  4.167
cat("Median Simpson: ", round(median(simpson_vals), 4), "\n")
#> Median Simpson:  0.0272
```

#### Model comparison via AIC and BIC

Compare nested models directly using
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html).

``` r
aic_vals <- c(
  Null     = AIC(fit_null),
  EnvHet   = AIC(fit_M1),
  `+Habitat` = AIC(fit_M2),
  `+Stream`  = AIC(fit_M3)
)
bic_vals <- c(
  Null     = BIC(fit_null),
  EnvHet   = BIC(fit_M1),
  `+Habitat` = BIC(fit_M2),
  `+Stream`  = BIC(fit_M3)
)
round(rbind(AIC = aic_vals, BIC = bic_vals), 2)
#>       Null EnvHet +Habitat +Stream
#> AIC 346.05 347.89   344.66  343.26
#> BIC 347.96 351.71   356.13  356.65
```

The model with the lowest AIC/BIC is preferred. Note that
`step.HubbellGLM()` can also be called directly for automated
backward/forward selection; see
[`?step.hubbell`](https://alessandrozito.github.io/HubbellGLM/reference/step.hubbell.md).

### Estimate $\sigma$

Before fitting any covariate model, estimate the accumulation-curve
growth rate from the null model.

``` r
sigma_hat <- estimate_sigma(formula = fit_M3$formula, data = BarroColorado)
#> Evaluating sigma =  0 
#> Evaluating sigma =  1.490116e-08 
#> Evaluating sigma =  -1 
#> Evaluating sigma =  -0.9999822 
#> Evaluating sigma =  -0.9063331 
#> Evaluating sigma =  -0.9063401 
#> Evaluating sigma =  -0.2471878 
#> Evaluating sigma =  -0.7269672 
#> Evaluating sigma =  -0.7269869 
#> Evaluating sigma =  -0.5742469 
#> Evaluating sigma =  -0.5742609 
#> Evaluating sigma =  -0.6462721 
#> Evaluating sigma =  -0.6462619 
#> Evaluating sigma =  -0.6373287 
#> Evaluating sigma =  -0.6373382 
#> Evaluating sigma =  -0.6373287
cat("Estimated sigma:", round(sigma_hat, 4), "\n")
#> Estimated sigma: -0.6373
```

A value near zero indicates near-logarithmic growth, which is consistent
with the Dirichlet-process baseline. Values in $(0,1)$ indicate
polynomial (heavier-tailed) growth, while $\sigma < 0$ implies that the
accumulation curve eventually converges to a finite number. —

## Session information

``` r
sessionInfo()
#> R version 4.3.1 (2023-06-16)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 20.04.6 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Berlin
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] HubbellGLM_1.0.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.36     desc_1.4.2        R6_2.5.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.0.6      knitr_1.51        htmltools_0.5.8.1
#>  [9] rmarkdown_2.27    lifecycle_1.0.4   cli_3.6.3         sass_0.4.9       
#> [13] pkgdown_2.2.0     textshaping_0.4.0 jquerylib_0.1.4   systemfonts_1.1.0
#> [17] compiler_4.3.1    rprojroot_2.0.3   rstudioapi_0.14   tools_4.3.1      
#> [21] ragg_1.5.2        bslib_0.8.0       evaluate_0.24.0   Rcpp_1.0.14      
#> [25] yaml_2.3.10       jsonlite_1.8.4    htmlwidgets_1.6.1 rlang_1.1.4      
#> [29] fs_1.6.0
```

------------------------------------------------------------------------

## References

Antoniak, Charles E. 1974. “Mixtures of Dirichlet Processes with
Applications to Bayesian Nonparametric Problems.” *The Annals of
Statistics* 2 (6): 1152–74.

Ewens, Warren J. 1972. “The Sampling Theory of Selectively Neutral
Alleles.” *Theoretical Population Biology* 3 (1): 87–112.

Ferguson, Thomas S. 1973. “A Bayesian Analysis of Some Nonparametric
Problems.” *The Annals of Statistics* 1 (2): 209–30.

Fisher, Ronald A., A. Steven Corbet, and C. B. Williams. 1943. “The
Relation Between the Number of Species and the Number of Individuals in
a Random Sample of an Animal Population.” *Journal of Animal Ecology* 12
(1): 42–58.

Hubbell, Stephen P. 2001. *The Unified Neutral Theory of Biodiversity
and Biogeography*. Princeton University Press.

Oksanen, Jari, Gavin L. Simpson, F. Guillaume Blanchet, Roeland Kindt,
Pierre Legendre, Peter R. Minchin, R. B. O’Hara, et al. 2024. *Vegan:
Community Ecology Package*. <https://CRAN.R-project.org/package=vegan>.

Pitman, Jim. 2006. *Combinatorial Stochastic Processes*. Vol. 1875.
Lecture Notes in Mathematics. Springer.

Rigon, Tommaso, Cheng-Long Hsu, and David B. Dunson. 2025. “A Bayesian
Theory for Estimation of Biodiversity.” *arXiv:2502.01333*.

Zito, Alessandro, Tommaso Rigon, Tomas Roslin, Pekka Niittynen, Paul D.
N. Hebert, Evgeny V. Zakharov, Sujeevan Ratnasingham, iBOL Consortium,
Otso Ovaskainen, and David B. Dunson. 2026. “Predicting Global
Biodiversity via Hubbell Regression.” *BioArxiv*.
