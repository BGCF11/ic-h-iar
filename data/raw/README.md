# Sentinel-2 input table

`Extracao_MataSantaGenebra_Integral_Corrigida.csv` contains the Phase B pixel/scene observations extracted through Google Earth Engine.

| Columns | Meaning |
| --- | --- |
| `longitude`, `latitude` | Pixel coordinates in decimal degrees |
| `B2`, `B3`, `B4`, `B8` | Surface reflectance on the `DN / 100` scale |
| `time` | Acquisition time in Unix milliseconds |

The table precedes same-day aggregation, chronological splitting and removal of deterministic components. Phase B performs those steps after loading the file.

See the [data dictionary](../../docs/DATA_DICTIONARY.md), [extraction settings](../metadata/README.md) and [execution commands](../../docs/REPRODUCING.md).
