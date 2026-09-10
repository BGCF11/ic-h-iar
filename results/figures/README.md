# Manuscript figures

| File | Content |
| --- | --- |
| `monte_carlo_max_bias.png` | Maximum componentwise absolute bias for each scenario and sample size |
| `monte_carlo_sd_scalar.png` | Monte Carlo standard deviation of the estimated scalar component `a` |
| `monte_carlo_time.png` | Mean estimation time by scenario and sample size |
| `map_temporal_persistence.png` | Estimated quaternion-norm surface |
| `map_vector_dominance.png` | Vector-dominance ratio |
| `map_b8_rmse.png` | Residual B8 RMSE |
| `map_persistence_gradient.png` | Persistence-surface Sobel gradient |

These are the manuscript display assets. The aggregate Monte Carlo values are in [results/tables](../tables/README.md), and the numerical georeferenced map values are in [results/geotiff](../geotiff/README.md).

The simulation script's `gerar_graficos_ic()` routine produces a separate three-panel diagnostic image with Portuguese labels. Both its bias and standard-deviation panels use scalar component `a`. The manuscript bias image above instead takes the maximum absolute bias across all four components. The diagnostic plot and the manuscript layout therefore represent different presentation steps.

The map script exports numerical rasters, while the manuscript PNGs also contain cartographic styling. It does not encode the final color ramps, clipping limits or page layout. PNG display images do not carry the numerical and spatial information of the GeoTIFFs.
