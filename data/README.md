# Data and results

The empirical input is a spatiotemporal table extracted from the [Harmonized Sentinel-2 SR collection](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED). Numerical result tables and georeferenced rasters are stored separately from the input and from the manuscript figure images.

## File catalogue

Paths below are relative to the repository root.

| File | Content |
| --- | --- |
| `data/raw/Extracao_MataSantaGenebra_Integral_Corrigida.csv` | Pixel/scene observations before local same-day aggregation, splitting and detrending |
| `results/tables/Resultados_Termografia_SantaGenebra_FaseB.csv` | Pixelwise quaternion estimates, training/test counts, predictive RMSEs and optimizer diagnostics |
| `results/tables/resultados_hiar_mc.csv` | Aggregate estimates, absolute biases, standard deviations, timings and iteration counts for the Monte Carlo experiment |
| `results/geotiff/Mapa_1_Resiliencia_Global.tiff` | Quaternion-norm surface: multispectral temporal persistence |
| `results/geotiff/Mapa_2_Dominancia_Vetorial_L1.tiff` | Ratio of the L1 magnitude of the quaternion vector part to the absolute scalar part |
| `results/geotiff/Mapa_3_RMSE_Infravermelho_B8.tiff` | One-step predictive RMSE in band B8 |
| `results/geotiff/Mapa_4_Gradiente_Morfologico.tiff` | Sobel-gradient magnitude of the persistence surface |

The original filenames and column names are used throughout the scripts. English definitions are given in the [data dictionary](../docs/DATA_DICTIONARY.md).

## Data access and execution

The seven numerical files are tracked with Git LFS. After cloning the repository with Git LFS installed, `git lfs pull` retrieves their contents. A small text file beginning with `version https://git-lfs.github.com/spec/v1` is a pointer, not a CSV or GeoTIFF dataset.

The research scripts resolve filenames from their execution directory. The [reproduction guide](../docs/REPRODUCING.md) copies the required input to a local directory under `runs/` before execution. Reference files remain in the catalogue locations above.

## Extraction and units

The extraction interval is `[2020-01-01, 2023-12-31)`. The mask retains `SCL == 4`, and the sampled bands are B2, B3, B4 and B8. Reflectance is stored as `DN / 100`, so a value of 20 represents physical reflectance 0.20. Acquisition time is in Unix milliseconds.

The input table precedes the chronological split and deterministic preprocessing. Multiple observations of a coordinate pair within a day are consolidated during Phase B, not during this export.

The [extraction metadata](metadata/README.md) and [study polygon](metadata/study_area.geojson) document the source settings. Re-extraction uses the live Earth Engine collection and can differ from the reference input if upstream data or processing change.

## Other script outputs

`h_iar_fase_a.py` produces `resultados_Fase_A_teste.csv` for the independent pilot. `fase_b_mapas.py` also produces `Mapa_5_Borda_Termodinamica_CSD.tiff`, a binary mask defined by quaternion norm at least 0.95. These are outputs of local execution, not additional inputs required by the workflow or entries in the seven-file reference set.

The Monte Carlo CSV contains aggregate statistics. Individual simulated series and per-repetition fits are not exported by the simulation script. Manuscript PNGs in [results/figures](../results/figures/README.md) are display assets; the tables and GeoTIFFs provide the numerical values.
