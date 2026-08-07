# Computational methods

## Analysis unit and input

Each row is one independent analysis unit, such as a guide-level aggregate.
The model input contains one numeric dose, one numeric response, and optionally
a categorical group. The code fails on missing, non-finite, or constant
responses by default. The Hill model additionally requires strictly positive
doses because it is linear in log dose.

Whether guide-level aggregates are statistically independent must be justified
from the experimental design. If biological replicates, batches, or repeated
measurements are present, these simple curves may not be the appropriate final
inferential model.

## Optional control normalization

For a control set \(C\), normalization follows the historical analysis:

\[
x_i^* = \frac{x_i}{\bar{x}_C}, \qquad
y_i^* = \frac{y_i}{\bar{y}_C}.
\]

Normalization is performed only when explicitly requested. The divisors are
written to `normalization_factors.csv`, and the transformed values are written
to `prepared_data.csv`.

## Polynomial candidates

For degree \(d\), `runmodelsaic_poly()` fits an ordinary least-squares model
using raw powers:

\[
y_i = \beta_0 + \sum_{j=1}^{d} \beta_j x_i^j + \epsilon_i,
\qquad \epsilon_i \sim N(0, \sigma^2).
\]

The only supported polynomial candidates are degrees 1 (linear) and 2
(quadratic). The function reports
coefficients, standard errors, t tests, R-squared, adjusted R-squared, Gaussian
log likelihood, AIC, AICc, BIC, RMSE, and the linear-versus-quadratic nested
F-test. Predictions remain within the observed dose range.

## Four-parameter log-logistic candidates

`runmodelsaic_drc()` retains the historical function name and fits the
log-logistic or Hill curve with `drc::drm()` and `drc::LL.4()`:

\[
f(x) = c + \frac{d-c}{1 + \exp\{b[\log(x)-\log(e)]\}},
\]

where \(b\) is the slope or Hill exponent, \(c\) is the lower asymptote,
\(d\) is the upper asymptote, and \(e\) is ED50. The sign of \(b\) determines
whether the curve increases or decreases.

The lower and upper asymptotes are fixed to the observed response minimum and
maximum of the fitted input rows (`drc::LL.4(fixed = c(NA, min(y), max(y), NA))`).
Only the slope and ED50 are estimated. Doses must be strictly positive. The
LL.4 likelihood, AIC, BIC, and free-coefficient standard errors are obtained
from the fitted `drc` object; AICc is calculated using its reported likelihood
degrees of freedom.

## Smoothing-spline candidates

`runmodelsaic_spline()` calls `stats::smooth.spline()` at a requested effective
degrees of freedom, with 3 as the default. The returned effective degrees of
freedom may differ slightly from the request. Spline AIC uses the effective
degrees of freedom for the mean curve plus one parameter for residual variance.
This is an approximation and should be interpreted together with curve
diagnostics and sensitivity analyses.

## Information criteria

For all candidate families, residuals are treated as independent Gaussian
errors. Given residual sum of squares RSS and \(n\) observations, the maximum
likelihood variance estimate is \(\widehat{\sigma}^2 = \mathrm{RSS}/n\). The
reported log likelihood is

\[
\ell = -\frac{n}{2}\left[\log(2\pi) + 1 +
\log(\widehat{\sigma}^2)\right].
\]

For total parameter count \(k\), including residual variance:

\[
\mathrm{AIC} = -2\ell + 2k,
\]

\[
\mathrm{AICc} = \mathrm{AIC} +
\frac{2k(k+1)}{n-k-1},
\]

when \(n > k+1\), and

\[
\mathrm{BIC} = -2\ell + k\log(n).
\]

AIC comparisons are valid only when candidates use the same response values,
observations, transformation, and likelihood definition. Delta AIC and Akaike
weights are descriptive model-selection summaries, not hypothesis tests or
evidence that a fitted mathematical shape is the biological mechanism.

## Reproducibility

The analysis runner records the exact prepared data, normalization factors,
command options, Git commit, dirty-worktree state, and R session information.
The example data generator fixes its random seed. The model-fitting workflow
uses no random numbers.
