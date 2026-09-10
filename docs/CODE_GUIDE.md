# Code guide

The scripts use Portuguese function and variable names. Comments and module/function docstrings are in English. This guide summarizes the interfaces using the exact identifiers expected by callers.

## Core: `src/python/h_iar.py`

| Function | Purpose |
| --- | --- |
| `calcular_matriz_transicao(phi_params, dt, F_tj)` | Fill a preallocated real 4 × 4 matrix representing left multiplication by the temporal quaternion power |
| `filtro_kalman_hiar(X, delta_t, phi_params, norma_phi_sq, R, X0, P0)` | Run Kalman prediction/update recursion; return innovations, inverse innovation covariances and their determinants |
| `calcular_nll_hiar(v, Lambda_inv_mat, det_Lambda_vec)` | Evaluate the Gaussian innovation negative log-likelihood, excluding the initial slot |
| `funcao_objetivo(phi_params_flat, X, delta_t, R, X0, P0)` | Evaluate the NLL using radial projection and an exterior quadratic penalty |
| `estimar_hiar(X, t, R_in=None)` | Centre data, set plug-in covariance scales, run L-BFGS-B and return fit diagnostics |

### Estimator inputs

- `X`: real numerical array of shape `(N, 4)`, with a fixed meaning and unit for each column. The empirical ordering is B2, B3, B4, B8.
- `t`: numerical array of shape `(N,)`, ordered chronologically in a consistent unit. The empirical application uses days. The routine computes `np.diff(t)` internally.
- `R_in`: observation-error covariance of shape `(4, 4)`. Default `None` selects `1e-6 * I4`; the empirical scripts explicitly supply `4 * I4`.
- `F_tj`: a reusable floating-point `(4, 4)` array, modified in place by the transition function.

Inputs must already be finite and appropriately prepared. The core has no general input-validation layer for arbitrary callers. Changing time units changes the interpretation of the estimated quaternion power.

### Estimator outputs

| Dictionary key | Meaning |
| --- | --- |
| `phi_otimo` | Final four-component parameter, radially projected if necessary |
| `nll_minima` | Optimizer-reported objective value; can include the exterior penalty and precede final output projection |
| `Sigma_otimo` | Approximate unit-interval process-noise matrix, `P0 * (1 - norm(phi_otimo)**2)` |
| `Rho_xi` | Post-fit correlation matrix of noninitial innovations; not a fitted covariance parameter in the likelihood |
| `convergencia` | SciPy optimizer `success` flag |
| `mensagem` | Optimizer termination message |
| `iteracoes` | Optimizer iteration count |
| `nfev` | Objective-function evaluation count |

The returned dictionary does not contain a standardized prediction API or the full state trajectory. Phase B constructs its own one-step residual predictions using `calcular_matriz_transicao`.

## Monte Carlo: `h_iar_monte_carlo.py`

| Function | Purpose |
| --- | --- |
| `gerar_tempos_biar(N)` | Generate cumulative observation times and gaps from the two-exponential mixture |
| `gerar_serie_hiar(N, t, dt_array, phi_real)` | Generate an isotropic synthetic H-IAR series under the specified true quaternion |
| `executar_monte_carlo_hiar()` | Execute the full simulation grid and return an aggregate pandas table |
| `gerar_graficos_ic(df_resultados)` | Plot scalar-component bias and SD, and mean execution time |

The script calls these functions at module scope. **Importing it launches the experiment.**

## Empirical scripts

| Script / function | Purpose |
| --- | --- |
| `fase_b_dados.py` / `exportar_malha_florestal_to_drive()` | Prepare the full study-area Sentinel-2 table and submit a Drive export |
| `h_iar_fase_a.py` / `extrair_dados_gee()` | Extract a small pilot window through synchronous `getRegion()` |
| `h_iar_fase_a.py` / `desazonalizar_e_destrendizar(df_agregado)` | Remove the fitted annual harmonic and trend from the pilot series |
| Both empirical fitting scripts / `modelo_fourier_tendencia(...)` | Evaluate intercept + annual cosine + annual sine + linear trend |
| `h_iar_fase_b.py` / `obter_parametros_sazonais(t, y)` | Fit the deterministic model to one training band; return `None` on `RuntimeError` |
| `h_iar_fase_b.py` / `processar_pixel(args)` | Process `(longitude, latitude, pixel_dataframe)` and return a result dictionary or `None` |
| `fase_b_mapas.py` | Rasterize the Phase B table and write the five GeoTIFF outputs |

The nested `processar_imagem` functions implement masking, band scaling and metadata preparation in their respective Earth Engine workflows. Consult the sources for the full signatures and the [data dictionary](DATA_DICTIONARY.md) for their tabular interfaces.

## Common Portuguese terms

| Preserved identifier fragment | English meaning |
| --- | --- |
| `treino` / `teste` | training / test |
| `tempo_dias` / `delta_t` | time in days / observation gap |
| `banda` / `bandas` | band / bands |
| `variancias_marginais` | marginal variances |
| `phi_otimo` / `phi_real` | estimated / true quaternion parameter |
| `norma_phi` | quaternion Euclidean norm |
| `erros_quadraticos` | squared errors |
| `resultados` / `iteracoes` | results / iterations |
| `resiliencia` | historical persistence-related label; see the methodological qualifications |

Legacy terms inside identifiers, filenames, CSV columns and plot strings are not independent ecological claims. Use the more precise definitions in [METHOD_NOTES.md](METHOD_NOTES.md).
