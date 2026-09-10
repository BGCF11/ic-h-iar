# Reproduction guide

## Environment and data access

Obtain the repository and its numerical files with Git LFS as described in the [README](../README.md). Create a Python virtual environment and install `requirements.txt`. The optional Earth Engine drivers additionally require `requirements-earthengine.txt` and authorized service access.

The dependency files list direct requirements rather than a historical environment lock. Original Python, package and BLAS/LAPACK versions are not recorded. Numba compiles numerical kernels on first use; optimizer behavior and elapsed time can depend on library versions, hardware and compilation overhead.

The commands below are run from the repository root unless stated otherwise. They work in Bash/Git Bash; `mkdir`, `cd` and `cp` also have corresponding aliases in PowerShell. Use a fresh local run directory for each execution.

## Script behavior

| Script | Execution behavior | File interface |
| --- | --- | --- |
| `h_iar.py` | Defines the numerical engine | No research-data file I/O |
| `h_iar_monte_carlo.py` | Starts the complete experiment, including on import | Writes aggregate CSV and diagnostic PNG |
| `fase_b_dados.py` | Authenticates and submits an Earth Engine export, including on import | Writes CSV to Google Drive |
| `h_iar_fase_a.py` | Authenticates, extracts and fits the pilot, including on import | Writes the Phase A CSV |
| `h_iar_fase_b.py` | Main execution is protected by `if __name__ == '__main__'` | Reads input CSV and writes pixelwise CSV |
| `fase_b_mapas.py` | Reads results and exports rasters, including on import | Reads pixelwise CSV and writes five GeoTIFFs |

Only the core is intended as an estimator interface. Importing all drivers as a test launches computations or service requests. Drivers resolve paths relative to the working directory and overwrite their fixed output filenames on repeated execution.

## 1. Monte Carlo experiment

```bash
mkdir runs/monte_carlo_local
cd runs/monte_carlo_local
python ../../src/python/h_iar_monte_carlo.py
cd ../..
```

The design uses four true parameter vectors, sample sizes 30, 100 and 300, and 1,000 repetitions per combination: 12,000 fits. The script sets `np.random.seed(42)` and uses NumPy random generation in ordinary Python functions.

| Local output | Content |
| --- | --- |
| `runs/monte_carlo_local/resultados_hiar_mc.csv` | 48 aggregate rows when all combinations complete |
| `runs/monte_carlo_local/graficos_hiar_mc.png` | Three-panel diagnostic plot: scalar bias, scalar SD and mean fit time |

The reference table is `results/tables/resultados_hiar_mc.csv`. Compare rows by `Caso`, `N` and `Parâmetro`. The [data dictionary](DATA_DICTIONARY.md) defines the statistics, including the `ddof=0` convention for SD. The script exports aggregates, not individual simulated series or fit records.

The fixed seed identifies the random experiment; it does not make optimization results or timing independent of the software environment. Timing covers estimation calls and can include first-use compilation. The manuscript's maximum-component bias image is described separately in [figure documentation](../results/figures/README.md).

## 2. Regional Phase B estimation and raster export

```bash
mkdir runs/phase_b_local
cp data/raw/Extracao_MataSantaGenebra_Integral_Corrigida.csv runs/phase_b_local/
cd runs/phase_b_local
python ../../src/python/h_iar_fase_b.py
python ../../src/python/fase_b_mapas.py
cd ../..
```

These steps run locally without Earth Engine. The first driver:

1. groups observations by coordinate pair and averages same-day records;
2. orders the series and assigns the final 10% to testing;
3. retains pixels with at least 30 training observations;
4. fits annual seasonality and linear trend to the training segment only;
5. estimates H-IAR on training residuals;
6. computes one-step test RMSE from the preceding observed residual.

The result is `runs/phase_b_local/Resultados_Termografia_SantaGenebra_FaseB.csv`. Processing uses `ProcessPoolExecutor()` with its default worker count and holds pixel groups and task objects in memory. Resource requirements depend on the input table and execution environment.

The second driver reads that CSV and writes Maps 1–5 in the same run directory. Maps 1–4 correspond to the reference rasters under `results/geotiff/`. Map 5 is an additional exploratory threshold mask.

Compare pixelwise results by `Longitude` and `Latitude`, since asynchronous processing can change row order. Numerical-success flags remain in the output, and the map script does not exclude rows solely because `Convergencia=False`. Worker exceptions return `None` and are omitted without a separate failure log.

## 3. Raster export from the reference table

This path generates rasters without repeating the parameter fits:

```bash
mkdir runs/maps_local
cp results/tables/Resultados_Termografia_SantaGenebra_FaseB.csv runs/maps_local/
cd runs/maps_local
python ../../src/python/fase_b_mapas.py
cd ../..
```

The local outputs are distinct files from the reference GeoTIFFs. Raster comparisons concern pixel values, NoData masks, coordinate system and transform; the map script does not implement manuscript cartographic styling.

## 4. Earth Engine extraction

This is an alternative route to obtain observations from the live source collection. It is not required for the local workflow above.

```bash
python -m pip install -r requirements-earthengine.txt
```

Both Earth Engine drivers call `ee.Authenticate()` and initialize project `hiar-497604`. Unmodified execution requires authorized access to that project. For another authorized project, a local adaptation of the initialization is required. The reference input CSV provides a service-independent entry point for regional analysis.

With project access configured:

```bash
mkdir runs/extraction_local
cd runs/extraction_local
python ../../src/python/fase_b_dados.py
cd ../..
```

The export is asynchronous and targets Google Drive folder `IC_Sensoriamento_H_IAR`, with description `Extracao_MataSantaGenebra_Integral_Corrigida`. After the task completes, download the CSV into a local run directory. Submission of the task alone does not imply completion of the export.

The code selects `[2020-01-01, 2023-12-31)`, with an [exclusive end date](https://developers.google.com/earth-engine/apidocs/ee-imagecollection-filterdate). Band values are already in percentage-reflectance units and `time` is in Unix milliseconds. Settings are documented in [extraction metadata](../data/metadata/README.md). See also [Earth Engine authentication](https://developers.google.com/earth-engine/guides/auth).

## 5. Independent Phase A pilot

With the same Earth Engine authorization:

```bash
mkdir runs/phase_a_local
cd runs/phase_a_local
python ../../src/python/h_iar_fase_a.py
cd ../..
```

The driver samples a 60 m buffer and processes at most 40 pixels, writing `resultados_Fase_A_teste.csv` in the run directory. It fits the full pilot series, unlike the Phase B holdout evaluation. Phase A is not a dependency of the regional workflow.

## Comparing executions

For a local execution, record the repository commit, input file hashes, Python and package versions, hardware and process/thread settings. A dependency snapshot can be saved with `python -m pip freeze`; it describes that execution's environment.

Numerical comparison uses the reference CSVs and GeoTIFFs. Exact reconstruction of the manuscript display layout is a separate presentation task: the diagnostic plotting routine and raster exporter do not encode every final layout or GIS styling choice.
