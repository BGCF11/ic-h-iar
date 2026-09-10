# Georeferenced maps

| File | Quantity |
| --- | --- |
| `Mapa_1_Resiliencia_Global.tiff` | Estimated quaternion norm |
| `Mapa_2_Dominancia_Vetorial_L1.tiff` | Vector-dominance ratio |
| `Mapa_3_RMSE_Infravermelho_B8.tiff` | Predictive RMSE in B8, in percentage points |
| `Mapa_4_Gradiente_Morfologico.tiff` | Spatial Sobel-gradient magnitude of the persistence surface |

The map script writes single-band `float32` GeoTIFFs in EPSG:4326, using NaN as NoData. Its angular cell size is 0.00009 degrees, with a half-cell correction of the upper-left origin. This corresponds to an approximate 10 m geographic grid; metric cell dimensions vary with latitude.

The first three quantities are assigned directly to raster cells. The gradient uses temporary nearest-neighbour filling followed by restoration of the valid-data mask. See the [data dictionary](../../docs/DATA_DICTIONARY.md) and [method notes](../../docs/METHOD_NOTES.md).

A local run of `fase_b_mapas.py` also produces `Mapa_5_Borda_Termodinamica_CSD.tiff`, an additional exploratory high-persistence mask. The four files above are the numerical rasters corresponding to the manuscript maps. Display images and their interpretation are documented in [figures](../figures/README.md).
