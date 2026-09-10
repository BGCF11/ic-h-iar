# Method notes

These notes describe the numerical specification, predictive evaluation and interpretation of the H-IAR implementation. Historical Portuguese names in the interfaces use the definitions given here.

## Covariances and quaternion dynamics

The core sets `P0 = diag(np.var(X, axis=0))` after componentwise centring. A diagonal covariance is not necessarily isotropic: isotropy additionally requires equal diagonal entries. The general transition is the real matrix of **left** quaternion multiplication and belongs to a restricted scaled-rotation family, not an unrestricted 4 × 4 vector autoregression.

The implemented state-noise covariance is `Q_j = P0 * (1 - ||Phi||**(2*dt))`. For admissible norms and nonnegative empirical variances it is diagonal nonnegative. It uses empirical marginal variances as scales; it does not generally preserve an anisotropic `P0` as an exact stationary covariance. Equality with `P0 - F_j @ P0 @ F_j.T` requires the relevant invariance condition, including the isotropic case. Monte Carlo **generation** uses `P0 = I4`; each fit still recomputes empirical diagonal scales.

The observation covariance is fixed to `4 * I4` in the empirical scripts: an assumed standard deviation of two percentage points per band. It is not inferred from those observations, and no sensitivity analysis is included. Monte Carlo fitting instead supplies `1e-6 * I4`; the generator adds state noise but no separate observation noise.

## Nearly real negative parameters

When `sqrt(b*b + c*c + d*d) < 1e-12`, the transition implementation sets the scalar power to `a**dt` for positive `a`, and to zero otherwise; its vector coefficients are zero. The nonpositive branch is an implementation convention, not a general definition of a noninteger power of a negative real quaternion.

The optimizer trajectory is not exported, so the effect of this branch cannot be assessed from aggregate outputs alone. The four Monte Carlo truth vectors have nonzero vector parts.

## Numerical termination and retained outputs

L-BFGS-B uses componentwise bounds `[-0.99, 0.99]`, `maxiter=2000` and `ftol=1e-9`; other optimizer settings are library defaults. Objective evaluation projects parameters outside `||Phi||^2 <= 0.99` and adds a quadratic penalty. Returned parameters are projected again when necessary.

`convergencia` / `Convergencia` is SciPy's `success` flag. It is not a guarantee of a global optimum. The Phase B script retains non-`None` worker returns irrespective of that flag. The map script drops missing parameter norms but does not filter by convergence.

A failed success flag differs from a pixel-processing exception: the former can still yield a retained parameter vector, whereas the latter returns `None` and is omitted. The supplied code does not write a separate log of omitted pixels. Monte Carlo summaries likewise do not exclude fits based on the success flag or save it per repetition.

## Predictive evaluation

Phase B fits deterministic seasonality/trend and the quaternion parameter using the training segment. Each held-out prediction is `F_dt @ previous_observed_residual`. It does not propagate the Kalman-filtered latent state and is not a recursively accumulated multi-step forecast. Previous test observations become available for the next one-step evaluation; they are not used to refit the parameters.

RMSE is computed against the next residual, in percentage points. Adding the same estimated deterministic component back to the prediction and target preserves that difference. A median of pixelwise RMSEs and a pooled error over all observations are different summary statistics.

## Ecological and spatial interpretation

`Resiliencia_Global` is retained as an interface name for the quaternion norm. Its operational interpretation is multispectral temporal persistence at a one-day reference interval. It is not a direct physiological resilience measurement. One parameter is estimated per pixel from the training segment; this alone does not establish a temporal early-warning trend.

Vector dominance is an exploratory, basis-dependent ratio. The `1e-8` denominator offset prevents literal division by zero but does not prevent very large values near a zero scalar component. It is not a percentage or a calibrated measure of causal coupling.

The Sobel gradient highlights spatial contrasts in the persistence surface. Temporary nearest-neighbour filling enables convolution, and the original mask is restored afterwards. Gradients within the mask can still depend on filled neighbours. These surfaces identify candidate boundaries, not independently validated ecological classes. No supervised classification, field-validation dataset or spatial covariance model is included.

The map code assigns coordinates to a grid using an angular spacing of 0.00009 degrees. It does not perform final-value interpolation for the first three maps. Multiple coordinate pairs can map to the same array cell; direct assignment does not aggregate such collisions. The Sobel response is in array-grid units and is not normalized by distance in metres.

## Monte Carlo definitions and computational reproducibility

The aggregate `SD` uses NumPy's default `ddof=0` (denominator `M`), and `Abs Bias` is the absolute bias of the mean estimate. The diagnostic script plots **scalar-component** bias and SD. It does not implement a maximum-across-components bias plot, even though such a panel is present among the English manuscript assets.

Exact original dependency versions were not retained. New runs may differ because of numerical libraries, optimizer defaults, hardware, process scheduling and compilation overhead. For numerical comparison, use the reference output tables and identify the source commit and local execution environment.
