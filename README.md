# H-IAR: quaternion autoregression for irregular multispectral time series

Research code, data and results for the **Hypercomplex Irregular Autoregressive (H-IAR)** model, with a Sentinel-2 application to the Mata de Santa Genebra ARIE, Brazil.

H-IAR represents four contemporaneous measurements as a quaternion and models their evolution over irregular observation intervals. The parameter is estimated through a Kalman-filter innovation likelihood. The empirical application uses surface reflectance in bands B2, B3, B4 and B8, in that order.

## Documentation

- [Reproduction guide](docs/REPRODUCING.md): installation, execution order, working directories and outputs.
- [Code guide](docs/CODE_GUIDE.md): estimator interface and script responsibilities.
- [Data and results](data/README.md): input table, numerical outputs and georeferenced maps.
- [Data dictionary](docs/DATA_DICTIONARY.md): original column names, units and derived quantities.
- [Method notes](docs/METHOD_NOTES.md): covariance specification, predictive evaluation and interpretation.

## Repository contents

| Location | Contents |
| --- | --- |
| `src/python/` | H-IAR estimator and five experiment, extraction, processing and mapping scripts |
| `data/raw/` | Sentinel-2 spatiotemporal input table |
| `data/metadata/` | Study polygon and extraction settings |
| `results/tables/` | Pixelwise Phase B results and aggregate Monte Carlo results |
| `results/geotiff/` | Four georeferenced maps used in the manuscript |
| `results/figures/` | Seven manuscript figure images |
| `docs/` | Technical documentation and methodological bibliography |
| `runs/` | Git-ignored working directories for local executions |

The research scripts have English comments and module/function docstrings. Identifiers, numeric constants, file interfaces and other string literals retain their original form.

## Obtain the code and data

The three CSVs and four GeoTIFFs use [Git Large File Storage](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage). With Git and Git LFS installed:

```bash
git lfs install
git clone https://github.com/BGCF11/ic-h-iar.git
cd ic-h-iar
git lfs pull
```

Use this workflow to obtain the full numerical files. A GitHub-generated ZIP can contain LFS pointers instead of the data, depending on the repository's archive settings.

## Install dependencies

From the repository root:

```bash
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` in Windows PowerShell, then run:

```bash
python -m pip install -r requirements.txt
```

Earth Engine extraction additionally uses `requirements-earthengine.txt`. Local Monte Carlo, Phase B fitting from the input CSV, and raster export do not require an Earth Engine account.

The dependency files list direct requirements. Exact historical environment versions are not recorded; numerical results and execution times can vary between environments. Execution commands and environment details are in the [reproduction guide](docs/REPRODUCING.md).

## Workflow

| Task | Script | Input | Output |
| --- | --- | --- | --- |
| Monte Carlo experiment | `h_iar_monte_carlo.py` | Four parameter scenarios; seed 42 | Aggregate CSV and diagnostic plot |
| Sentinel-2 extraction | `fase_b_dados.py` | Earth Engine collection and study polygon | Spatiotemporal CSV exported to Google Drive |
| Regional estimation | `h_iar_fase_b.py` | Spatiotemporal CSV | Pixelwise parameter estimates and test RMSEs |
| Raster export | `fase_b_mapas.py` | Pixelwise result CSV | Five GeoTIFFs, including one additional exploratory mask |
| Pilot analysis | `h_iar_fase_a.py` | Small Earth Engine sampling window | Phase A result CSV |

Phase A is an independent pilot. Phase B can run directly from the input table in `data/raw/`. The four manuscript rasters are in `results/geotiff/`; the fifth raster is an additional output of the map script.

The drivers use filenames relative to the **current working directory**. Follow the commands in the reproduction guide to run them under `runs/` and keep the reference data and results separate. Several drivers execute work when imported.

## Scientific interpretation

Phase B uses a chronological 90/10 training/test split. Annual seasonality and linear trend are fitted to the training segment, and H-IAR is estimated on its residuals. Test RMSE measures one-step propagation from the preceding observed residual.

`Resiliencia_Global` is the historical column name for the estimated quaternion norm, interpreted as a temporal-persistence descriptor. The vector-dominance ratio and spatial Sobel gradient are exploratory descriptors. Their definitions and the distinction between numerical convergence and ecological interpretation are given in the [method notes](docs/METHOD_NOTES.md).

## Citation and research team

Citation metadata are provided in [CITATION.cff](CITATION.cff). A citation should identify the commit or archived release used.

| Author | Affiliation |
| --- | --- |
| Bruno Gonçalves C. Filho | School of Mechanical Engineering (FEM), University of Campinas (UNICAMP), Campinas, Brazil |
| Angelo Calil Bianchi | Institute of Science and Technology (ICT), Federal University of São Paulo (UNIFESP), São José dos Campos, Brazil |
| Aluísio de Souza Pinheiro | Department of Statistics, Institute of Mathematics, Statistics and Scientific Computing (IMECC), UNICAMP, Campinas, Brazil |

Supported by FAPESP grants **2023/02538-0** and **2025/21329-8**. Regional processing used the Zurich computing resources at IMECC/UNICAMP.

Source observations: [Copernicus Sentinel-2 SR Harmonized, Google Earth Engine catalogue](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED). Methodological references are collected in [docs/references.bib](docs/references.bib).
