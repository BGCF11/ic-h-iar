"""
Phase A — Pointwise empirical H-IAR calibration with Sentinel-2 data
=================================================================
This script implements Phase A of the H-IAR empirical application in the
Mata de Santa Genebra ARIE. It extracts multispectral Sentinel-2 Level-2A
Surface Reflectance series through Google Earth Engine, retains pixels
classified as vegetation by the SCL band, removes deterministic seasonal
and trend components, and passes the residuals to the H-IAR estimator.

Methodological purpose of Phase A:
evaluate the physical-statistical workflow on a controlled pilot area
before full regional processing in Phase B.

References:
- Google Earth Engine Data Catalog: COPERNICUS/S2_SR_HARMONIZED.
- Google Earth Engine API: filterBounds, filterDate, updateMask, getRegion.
- Wilson, Knight & McRoberts (2018): harmonic regression in remote-sensing
  time series.
- SciPy: curve_fit for nonlinear least-squares fitting.
- Elorrieta et al. (2019, 2021): CIAR/BIAR and irregular autoregressive estimation.
"""

import ee
import pandas as pd
import numpy as np
from h_iar import estimar_hiar

# Authenticate and initialize Google Earth Engine for the research project.
ee.Authenticate()
ee.Initialize(project='hiar-497604')

def extrair_dados_gee():
    """
    Extract multispectral Sentinel-2 Level-2A Surface Reflectance time series
    from the harmonized COPERNICUS/S2_SR_HARMONIZED collection for a controlled
    spatial window in the Mata de Santa Genebra.

    Processing steps:
    1. retain pixels classified as vegetation by SCL == 4;
    2. select bands B2, B3, B4 and B8;
    3. rescale digital numbers to percentage reflectance;
    4. preserve the system:time_start timestamp;
    5. return a spatiotemporal table for the H-IAR workflow.

    Reference: Google Earth Engine Data Catalog — Sentinel-2 SR Harmonized;
    SCL class table.
    """
    # Centre of the Phase A sampling window in Mata de Santa Genebra.
    ponto_central = ee.Geometry.Point([-47.115264, -22.819975])

    # Circular extraction window with a 60 m radius around the centre.
    # Phase A uses a small area for pilot evaluation of the pipeline.
    area_estudo = ponto_central.buffer(60)

    def processar_imagem(imagem):
        # Select the Scene Classification Layer (SCL) from Sentinel-2 L2A.
        # SCL class 4 corresponds to "Vegetation" in the official catalogue.
        scl = imagem.select('SCL')

        # Retain only pixels classified as vegetation.
        mascara_vegetacao = scl.eq(4)
        
        # Select the 10 m optical bands used as H-IAR vector components:
        # B2 (blue), B3 (green), B4 (red) and B8 (near infrared).
        #
        # In GEE, Sentinel-2 SR bands have a scale factor of 0.0001;
        # physical reflectance = DN / 10000.
        # Dividing by 100 converts DN to percentage reflectance:
        # DN/100 = 100 * physical reflectance, approximately on a [0, 100] scale.
        #
        # Percentage scaling prevents 4x4 covariance determinants from being
        # artificially close to the H-IAR diagonal-jitter threshold solely due to units.
        bandas = imagem.select(['B2', 'B3', 'B4', 'B8']).divide(100)
        
        # Apply the vegetation mask and retain system:time_start,
        # needed to reconstruct each pixel's irregular time series.
        return bandas.updateMask(mascara_vegetacao).copyProperties(imagem, ['system:time_start'])

    # Harmonized Sentinel-2 Level-2A Surface Reflectance collection.
    # HARMONIZED corrects the radiometric offset associated with
    # PROCESSING_BASELINE >= 04.00, aligning newer and older scenes.
    #
    # filterBounds restricts the collection to the Phase A window;
    # filterDate selects [2020-01-01, 2023-12-31), with an exclusive end date;
    # map applies the SCL mask and percentage scaling to each image.
    colecao = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
               .filterBounds(area_estudo)
               .filterDate('2020-01-01', '2023-12-31')
               .map(processar_imagem))

    # Spatiotemporal extraction at 10 m scale.
    # getRegion returns rows containing id, longitude, latitude, time and band
    # values for each pixel-image tuple within the region.
    extraidos = colecao.getRegion(area_estudo, 10).getInfo()
    
    # Convert the GEE response to a table, separating the header from observations.
    cabecalho = extraidos[0]
    dados = extraidos[1:]
    df = pd.DataFrame(dados, columns=cabecalho)
    
    # Remove observations missing any of the four H-IAR bands.
    # These NaNs arise from masking and observational gaps.
    df = df.dropna(subset=['B2', 'B3', 'B4', 'B8']).reset_index(drop=True)
    
    return df

# Extract the pilot sampling window.
df_bruto = extrair_dados_gee()

from scipy.optimize import curve_fit

def modelo_fourier_tendencia(t, beta_0, beta_1, beta_2, beta_3):
    """
    Evaluate the deterministic model of annual seasonality and linear trend.

    The annual harmonic component is:
    beta_1*cos(omega*t) + beta_2*sin(omega*t), with omega = 2*pi/365.25.

    The beta_3*t term represents long-term linear drift. Removing the fitted
    trend before H-IAR estimation separates this deterministic component from
    the dynamics described by the autoregressive parameter Phi.

    References: Fourier harmonic regression in remote-sensing time series;
    SciPy curve_fit for nonlinear least-squares estimation.
    """
    omega = 2 * np.pi / 365.25
    return beta_0 + beta_1 * np.cos(omega * t) + beta_2 * np.sin(omega * t) + beta_3 * t

def desazonalizar_e_destrendizar(df_agregado):
    """
    Remove the fitted deterministic components and return the residuals.
    """
    df_limpo = df_agregado.copy()
    t = df_limpo['tempo_dias'].values
    
    for banda in ['B2', 'B3', 'B4', 'B8']:
        y = df_limpo[banda].values
        
        # Initial values for the deterministic fit:
        # beta_0 is the band mean; harmonic and trend coefficients start at zero.
        p0 = [np.mean(y), 0.0, 0.0, 0.0] 
        
        parametros_otimos, _ = curve_fit(modelo_fourier_tendencia, t, y, p0=p0)
        y_teorico = modelo_fourier_tendencia(t, *parametros_otimos)
        
        # Remove annual seasonality and the linear trend,
        # retaining the irregular residual component supplied to H-IAR.
        df_limpo[banda] = y - y_teorico
        
    return df_limpo

# Common time reference: 2020-01-01 in days since the Unix epoch.
# Use a common time origin before computing Delta_t.
TEMPO_ZERO_DIAS = 1577836800000 / (1000 * 60 * 60 * 24)

resultados_fase_A = []
pixeis_processados = 0

# Reconstruct each pixel series from its longitude/latitude pair.
for (lon, lat), df_pixel in df_bruto.groupby(['longitude', 'latitude']):
    
    # Phase A operational limit: process at most 40 pixels in this supplied script
    # to evaluate the pipeline before regional Phase B processing.
    if pixeis_processados >= 40:
        break
        
    df_pixel = df_pixel.copy()
    
    # 1. Convert system:time_start from milliseconds to days since 2020-01-01.
    df_pixel['tempo_dias'] = (df_pixel['time'] / (1000 * 60 * 60 * 24)) - TEMPO_ZERO_DIAS
    
    # 2. Merge within-day observations associated with orbital overlaps.
    # Average observations of the same pixel on the same day to avoid
    # near-zero Delta_t and redundant radiometry in the H-IAR filter.
    df_pixel['dia_inteiro'] = np.floor(df_pixel['tempo_dias'])
    df_agregado = df_pixel.groupby('dia_inteiro').mean(numeric_only=True).reset_index()
    
    # 3. Sort chronologically before H-IAR computes Delta_t internally.
    df_agregado = df_agregado.sort_values('tempo_dias')
    
    # Minimum temporal-information requirement:
    # discard pixels with fewer than 30 observations after masking and aggregation,
    # for example where cloud masking leaves insufficient observations.
    if len(df_agregado) < 30:
        continue
    
    # 4. Remove deterministic components before H-IAR estimation.
    # Removing seasonality and trend targets the stationary residual formulation;
    # this preprocessing alone does not establish stationarity.
    try:
        df_processado = desazonalizar_e_destrendizar(df_agregado)
    except RuntimeError:
        # Fallback if the harmonic fit raises a numerical RuntimeError:
        # remove only the band mean to retain a centred series.
        df_processado = df_agregado.copy()
        for b in ['B2', 'B3', 'B4', 'B8']:
            df_processado[b] = df_processado[b] - df_processado[b].mean()
            
    # 5. Assemble H-IAR inputs:
    # X_real contains the four spectral components of the quaternion vector;
    # t_real contains the corresponding irregular observation times.
    X_real = df_processado[['B2', 'B3', 'B4', 'B8']].values
    t_real = df_processado['tempo_dias'].values
    
    # A priori observation covariance on the percentage-reflectance scale.
    # An assumed standard deviation of 2 percentage points
    # is represented by variance 2.0^2 for each component.
    R_instrumental = np.eye(4) * (2.0**2)
    
    # 6. Estimate the H-IAR quaternion parameter Phi for the current pixel.
    res = estimar_hiar(X_real, t_real, R_in=R_instrumental)
    
    # Package pixelwise results.
    # Phi_Otimizado stores the four estimated components of Phi = (a,b,c,d).
    # Resiliencia_Global is the historical field name for ||Phi||_2 (persistence).
    resultados_fase_A.append({
    'ID_Pixel': pixeis_processados + 1,
    'Longitude': lon, 
    'Latitude': lat,
    'N_observacoes': len(t_real),
    'Phi_Otimizado': res['phi_otimo'].tolist(),
    'Resiliencia_Global': np.linalg.norm(res['phi_otimo']),
    'Convergencia': res['convergencia'],
    'Iteracoes': res.get('iteracoes', 'N/A') 
})
    
    pixeis_processados += 1

    df_resultados_A = pd.DataFrame(resultados_fase_A)

# Export the Phase A table for subsequent analysis and spatial comparison.
df_resultados_A.to_csv('resultados_Fase_A_teste.csv', index=False)
