"""
Phase B — Rasterization and GeoTIFF mapping of H-IAR results
==========================================================
This script converts the consolidated Phase B table into georeferenced
GeoTIFF maps of the Mata de Santa Genebra ARIE.

Rasterization maps coordinates directly to array indices, without spatial
interpolation of the final parameter and RMSE surfaces. Cells without
assigned observations remain missing. The gradient calculation uses
temporary nearest-neighbour filling and then restores the valid-data mask.

The derived quantities include:
- Vector dominance: the L1 magnitude of the vector part of Phi divided by
  the absolute scalar part, with a 1e-8 denominator offset;
- High-persistence mask: a binary threshold at ||Phi|| >= 0.95, retained
  under the historical CSD filename;
- Morphological gradient: Sobel magnitude of the persistence surface.

References:
- Rasterio: GeoTIFF writing and affine transformations.
- SciPy: griddata(method='nearest') and ndimage.sobel.
"""

import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import sobel
import rasterio
from rasterio.transform import from_origin

# ===================================================
# 1. READ THE H-IAR TABLE AND COMPUTE MAPPING METRICS
# ===================================================

print("Carregando resultados da Fase B...")

# Load the consolidated table produced by regional H-IAR processing.
# Each row represents a processed pixel with coordinates, Phi parameters,
# the parameter norm and bandwise predictive metrics.
df = pd.read_csv('Resultados_Termografia_SantaGenebra_FaseB.csv')

# Remove rows with a missing quaternion norm.
# Pixels without Resiliencia_Global cannot enter the derived raster maps.
df = df.dropna(subset=['Resiliencia_Global']).reset_index(drop=True)

# Vector dominance: an exploratory descriptor defined for this H-IAR study.
# Ratio of the L1 magnitude of the imaginary part of Phi to the absolute
# scalar component, with a numerical denominator offset:
#     (|b| + |c| + |d|) / (|a| + 1e-8)
# The 1e-8 term prevents division by zero when Phi_a is numerically zero.
massa_imaginaria = np.abs(df['Phi_b']) + np.abs(df['Phi_c']) + np.abs(df['Phi_d'])
df['Dominancia_Vetorial'] = massa_imaginaria / (np.abs(df['Phi_a']) + 1e-8)

# Geographic coordinates of the retained pixels.
# Used below for temporary nearest-neighbour filling
# to enable the Sobel convolution.
pontos = df[['Longitude', 'Latitude']].values

# =============================================
# 2. DIRECT RASTERIZATION WITHOUT FINAL-VALUE INTERPOLATION
# =============================================

print("Mapeando matrizes espaciais em resolução nativa...")


# Operational angular resolution of the output grid in the study area.
# Delta = 0.00009 degrees is approximately 10 m and is used to
# convert longitude/latitude into discrete array indices.
resolucao_graus = 0.00009

# Observed spatial extent of the retained grid.
# Bounds are computed from the pixels actually processed.
min_lon, max_lon = df['Longitude'].min(), df['Longitude'].max()
min_lat, max_lat = df['Latitude'].min(), df['Latitude'].max()

# Raster array dimensions.
# Discretize the continuous geographic extent into cells of size Delta.
# The +1 retains the extreme cells of the observed grid.
largura = int(np.round((max_lon - min_lon) / resolucao_graus)) + 1
altura = int(np.round((max_lat - min_lat) / resolucao_graus)) + 1

# Direct mapping from coordinates to array indices.
# indices_x increases west to east; indices_y increases north to south,
# following the usual raster orientation.
#
# No interpolation is used: each H-IAR value is assigned to the cell
# corresponding to its observed pixel coordinates.
indices_x = np.round((df['Longitude'] - min_lon) / resolucao_graus).astype(int)
indices_y = np.round((max_lat - df['Latitude']) / resolucao_graus).astype(int)

# Initialize with NaN.
# Cells without an assigned value remain missing, retaining the observed
# gaps and study-area outline without filling the final parameter surfaces.
matriz_resiliencia = np.full((altura, largura), np.nan)
matriz_dominancia = np.full((altura, largura), np.nan)
matriz_rmse_b8 = np.full((altura, largura), np.nan)

# Direct assignment of observed/estimated values.
# Each table row is assigned to one raster cell.
matriz_resiliencia[indices_y, indices_x] = df['Resiliencia_Global'].values
matriz_dominancia[indices_y, indices_x] = df['Dominancia_Vetorial'].values
matriz_rmse_b8[indices_y, indices_x] = df['RMSE_B8'].values

# ========================================================
# 3. DERIVED MAPS: HIGH-PERSISTENCE MASK AND MORPHOLOGICAL GRADIENT
# ========================================================

print("Processando tensores de gradiente...")

# A. Binary high-persistence mask; the legacy name refers to the CSD motivation.
# Pixels with ||Phi_hat|| >= 0.95 receive 1.0; the others receive 0.0.
# The 0.95 threshold is exploratory and does not independently establish CSD.
matriz_borda_csd = np.where(matriz_resiliencia >= 0.95, 1.0, 0.0)

# Preserve missing values outside the retained raster domain.
matriz_borda_csd[np.isnan(matriz_resiliencia)] = np.nan

# B. Morphological gradient of the temporal-persistence surface.
# Create a regular grid with the same extent and resolution as the raster.
grid_x, grid_y = np.meshgrid(
    np.arange(min_lon, max_lon + resolucao_graus/2, resolucao_graus),
    np.arange(max_lat, min_lat - resolucao_graus/2, -resolucao_graus)
)

# Match the grid dimensions to the raster arrays.
grid_x, grid_y = grid_x[:altura, :largura], grid_y[:altura, :largura]

# Temporary nearest-neighbour filling.
# Fill missing values temporarily rather than interpolating the final parameter map,
# allowing application of the Sobel convolution operator.
#
# Restore the original NaN mask after computing the gradient.
resiliencia_preenchida = griddata(pontos, df['Resiliencia_Global'].values, (grid_x, grid_y), method='nearest')


# Apply the Sobel operator along each array axis.
# axis=0 captures row-direction variation; axis=1 captures column-direction
# variation. The final magnitude is sqrt(dx^2 + dy^2).
dx = sobel(resiliencia_preenchida, axis=0)
dy = sobel(resiliencia_preenchida, axis=1)
matriz_gradiente = np.hypot(dx, dy)

# Restore the original missing-data mask.
# Remove gradient values outside the original valid domain after temporary
# filling; gradients within that domain can still depend on neighbouring fills.
matriz_gradiente[np.isnan(matriz_resiliencia)] = np.nan

# ==========================================
# 4. GEOREFERENCED GEOTIFF EXPORT
# ==========================================

print("Exportando cartografia em formato .TIFF...")

# Half-pixel correction of the cartographic origin.
# from_origin expects the upper-left raster corner; since the coordinates
# represent cell centres, shift the origin by Delta/2:
# minimum longitude - Delta/2 and maximum latitude + Delta/2.
origem_lon = min_lon - (resolucao_graus / 2)
origem_lat = max_lat + (resolucao_graus / 2)
transformacao = from_origin(origem_lon, origem_lat, resolucao_graus, resolucao_graus)

# Final collection of exported GeoTIFFs; historical filenames are retained.
# Map 1: global quaternion norm ||Phi_hat||.
# Map 2: L1 ratio of the imaginary part to the scalar component.
# Map 3: one-step-ahead B8 RMSE in percentage points.
# Map 4: morphological gradient of the temporal-persistence surface.
# Map 5: binary high-persistence mask defined by ||Phi_hat|| >= 0.95.
dicionario_mapas = {
    'Mapa_1_Resiliencia_Global.tiff': matriz_resiliencia,
    'Mapa_2_Dominancia_Vetorial_L1.tiff': matriz_dominancia,
    'Mapa_3_RMSE_Infravermelho_B8.tiff': matriz_rmse_b8,
    'Mapa_4_Gradiente_Morfologico.tiff': matriz_gradiente,
    'Mapa_5_Borda_Termodinamica_CSD.tiff': matriz_borda_csd
}

for nome_ficheiro, matriz in dicionario_mapas.items():
    # Convert to float32 before writing.
    # Reduce file size and improve GIS software compatibility,
    # with float32 precision for the cartographic metrics.
    matriz_float32 = matriz.astype(np.float32)
    
    # Write a georeferenced, single-band GeoTIFF.
    # crs='EPSG:4326' stores longitude/latitude in degrees;
    # transform specifies the origin and pixel size;
    # nodata=np.nan explicitly represents missing data.
    with rasterio.open(
        nome_ficheiro,
        'w',
        driver='GTiff',
        height=matriz_float32.shape[0],
        width=matriz_float32.shape[1],
        count=1,
        dtype=matriz_float32.dtype,
        crs='EPSG:4326',
        transform=transformacao,
        nodata=np.nan
    ) as ficheiro_tiff:
        ficheiro_tiff.write(matriz_float32, 1)

print("Operação concluída. Mapas termográficos gerados com sucesso.")
