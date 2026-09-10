# Data dictionary

Column spellings below are the exact strings used by the research scripts and numerical tables. Definitions and units are provided in English.

## Phase B input

File: `Extracao_MataSantaGenebra_Integral_Corrigida.csv`

Each row is a sampled pixel/scene observation after the Earth Engine vegetation mask, before local within-day aggregation. Repeated coordinate pairs and multiple records within a day are expected at this stage.

| Field | Meaning and unit |
| --- | --- |
| `longitude` | Pixel longitude, decimal degrees |
| `latitude` | Pixel latitude, decimal degrees |
| `B2` | Blue-band surface reflectance, percentage points of unit reflectance (DN / 100) |
| `B3` | Green-band surface reflectance, same scale |
| `B4` | Red-band surface reflectance, same scale |
| `B8` | Near-infrared surface reflectance, same scale |
| `time` | Scene `system:time_start`, Unix milliseconds |

The quaternion component order is `(B2, B3, B4, B8) → (scalar, i, j, k)`. These are not standardized z-scores or raw Sentinel-2 digital numbers. For example, a reflectance value of 20 on this scale denotes physical reflectance 0.20; an error of 2 is two percentage points, not a 2% relative error.

During Phase B, time is converted to days since 2020-01-01. Grouping uses the floor of that day count; the aggregated time is the mean acquisition time for the group. After sorting, `int(0.90 * n_total)` observations are used for training. At least 30 training observations are required; this is not simply a cutoff of 30 total dates.

## Phase B output

File: `Resultados_Termografia_SantaGenebra_FaseB.csv`

Each row is a pixel result returned by the worker; asynchronous completion means rows need not be spatially sorted.

| Field | Meaning and unit |
| --- | --- |
| `Longitude`, `Latitude` | Pixel coordinates, decimal degrees; note capitalization differs from the input |
| `N_Treino` | Number of training dates after within-day aggregation |
| `N_Teste` | Number of held-out test dates |
| `Phi_a`, `Phi_b`, `Phi_c`, `Phi_d` | Estimated scalar and three vector components of the quaternion parameter |
| `Resiliencia_Global` | Euclidean norm of the estimated quaternion; temporal-persistence factor at the one-day reference interval |
| `RMSE_B2`, `RMSE_B3`, `RMSE_B4`, `RMSE_B8` | Test one-step residual RMSE, in percentage points |
| `Convergencia` | Optimizer `success` Boolean; not a certificate of global optimality or scientific validation |
| `Iteracoes` | Number of L-BFGS-B iterations returned by the estimator |

Parameters are estimated on training residuals. Test predictions propagate the preceding observed residual (the last training residual for the first test point). The same fitted deterministic component can be added back to predictions without changing the paired errors, but the supplied script computes the errors in residual space.

`Convergencia=False` rows are not excluded by the map script. Missing residual RMSE values are initialized as NaN; pandas exports missing numerical fields as empty CSV fields by default. Pixel workers returning `None` have no row and no separately saved failure reason. The result table records the flags and missing fields for retained pixels.

## Monte Carlo aggregate output

File: `resultados_hiar_mc.csv`

| Field | Meaning and unit |
| --- | --- |
| `Caso` | Scenario name: `Caso 1 (+/+)`, `Caso 2 (-/-)`, `Caso 3 (-/+)`, `Caso 4 (+/-)` |
| `N` | Sample size: 30, 100 or 300 |
| `Parâmetro` | Quaternion component: `a`, `b`, `c` or `d` |
| `Valor Real` | True component value in the simulation |
| `Média Estimada` | Mean component estimate over 1,000 fits |
| `Abs Bias` | Absolute difference between the mean estimate and true value; not the mean absolute error of individual fits |
| `SD` | Standard deviation of component estimates computed by `np.std(..., axis=0)`, hence `ddof=0` and denominator 1,000 |
| `Tempo Total (s)` | Mean elapsed estimation time per fit, in seconds; despite the historical name, not the summed time of all repetitions |
| `Iterações` | Mean optimizer iteration count per fit |
| `Tempo/Iteração (ms)` | `1000 * mean_time / max(mean_iterations, 1)`, in milliseconds; not a mean of individual time/iteration ratios |

There are 48 aggregate rows if all combinations finish. Time and iteration summaries are repeated in each of the four component rows for a scenario/sample-size pair. The CSV does not store the 12,000 individual fits.

True parameter vectors in scenario order are `(0.7, 0.3, 0.3, 0.3)`, `(-0.7, -0.3, -0.3, -0.3)`, `(-0.9, 0.15, 0.15, 0.15)` and `(0.9, -0.15, -0.15, -0.15)`.

## Phase A local output

File: `resultados_Fase_A_teste.csv`

Fields are `ID_Pixel`, `Longitude`, `Latitude`, `N_observacoes`, `Phi_Otimizado`, `Resiliencia_Global`, `Convergencia`, `Iteracoes`. `Phi_Otimizado` is a four-element list serialized into one CSV field, not four separate columns. `N_observacoes` counts the full pilot series after same-day aggregation. The script processes at most 40 pixels per execution. This output is generated by the independent pilot and is not an input to Phase B.

## Derived raster quantities

Maps 1–4 are the reference rasters. Map 5 is an additional output generated by a local run of the mapping script.

| Raster | Stored quantity |
| --- | --- |
| Map 1 | `sqrt(Phi_a**2 + Phi_b**2 + Phi_c**2 + Phi_d**2)` |
| Map 2 | `(|Phi_b| + |Phi_c| + |Phi_d|) / (|Phi_a| + 1e-8)` |
| Map 3 | `RMSE_B8`, in percentage points |
| Map 4 | Magnitude of the two default SciPy Sobel responses on a temporarily nearest-filled persistence surface, followed by restoration of the valid-data mask |
| Map 5 | 1 if the parameter norm is at least 0.95, otherwise 0; NoData retained outside the valid domain |

Map 2 can become very large as the scalar component approaches zero; it is not a bounded proportion. Map 4 is an array-grid contrast measure, not a gradient normalized per metre. Maps use NaN NoData and the original historical Portuguese filenames.
