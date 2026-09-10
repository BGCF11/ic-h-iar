# Extraction metadata

| Setting | Value |
| --- | --- |
| Study area | Mata de Santa Genebra ARIE, Brazil |
| Collection | `COPERNICUS/S2_SR_HARMONIZED` |
| Start | `2020-01-01`, inclusive |
| End | `2023-12-31`, exclusive |
| Vegetation mask | Sentinel-2 `SCL == 4` |
| Band order | B2, B3, B4, B8 |
| Band scaling | Digital numbers divided by 100: percentage reflectance |
| Time field | `system:time_start`, exported as numerical `time` in Unix milliseconds |
| Phase B sampling | `scale=10`, `projection='EPSG:4326'`, `geometries=False`, `dropNulls=True` |
| Phase B export | CSV to Drive folder `IC_Sensoriamento_H_IAR` |
| Export description | `Extracao_MataSantaGenebra_Integral_Corrigida` |
| Earth Engine project in the scripts | `hiar-497604` |
| Phase A region | Point `[-47.115264, -22.819975]`, buffered by 60 m |
| Phase A extraction | `getRegion(area_estudo, 10).getInfo()` |

[study_area.geojson](study_area.geojson) contains the polygon vertices defined in `fase_b_dados.py`, with the ring explicitly closed. This is the operational study polygon, not a certified administrative boundary. The `ee.Geometry.Polygon` call in the extraction script defines the geometry used by Earth Engine.

The sampling scale and the output raster grid describe different steps. Earth Engine sampling uses a scale of 10 m; the map script later constructs a geographic grid with a spacing of 0.00009 degrees.

The two Earth Engine drivers authenticate and initialize the project shown above. Local analysis from the input CSV does not use that service. See the [reproduction guide](../../docs/REPRODUCING.md) for access requirements.

Sources: [Sentinel-2 SR Harmonized catalogue](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED) and [filterDate interval semantics](https://developers.google.com/earth-engine/apidocs/ee-imagecollection-filterdate).
