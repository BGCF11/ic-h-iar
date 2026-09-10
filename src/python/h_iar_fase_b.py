"""
Phase B — Regional pixelwise H-IAR processing and predictive evaluation
=====================================================================
This script processes the complete Mata de Santa Genebra spatiotemporal
dataset previously exported through Google Earth Engine. The table is
grouped by pixel, and each time series is processed independently in parallel.

For each pixel:
1. convert raw timestamps to days since 2020-01-01;
2. merge within-day observations by averaging;
3. apply a chronological 90/10 split before any deterministic fitting;
4. fit annual seasonality and linear trend using training data only;
5. estimate the quaternion parameter Phi with the H-IAR engine;
6. evaluate one-step-ahead test RMSE for B2, B3, B4 and B8.

References:
- Elorrieta et al. (2019, 2021): CIAR/BIAR and irregular autoregressive estimation.
- Kuipers (1999) and Mebius (2005): quaternion operations used by the H-IAR engine.
- Kapoor & Narayanan (2023): data leakage in predictive workflows.
- Hyndman & Athanasopoulos: temporal validation and forecast evaluation.
- SciPy: curve_fit for nonlinear least squares.
- Python concurrent.futures: ProcessPoolExecutor for parallel execution.
"""

import pandas as pd
import numpy as np
import time
import concurrent.futures
from scipy.optimize import curve_fit

# Import the H-IAR engine:
# estimar_hiar estimates Phi on the training data;
# calcular_matriz_transicao builds F_{Delta_t} for one-step-ahead evaluation.
from h_iar import estimar_hiar, calcular_matriz_transicao

def modelo_fourier_tendencia(t, beta_0, beta_1, beta_2, beta_3):
    """
    Evaluate the deterministic model for annual seasonality and linear drift.

    y(t) = beta_0 + beta_1*cos(omega*t) + beta_2*sin(omega*t) + beta_3*t,
    with omega = 2*pi/365.25.

    The harmonic component represents the annual vegetation cycle; the linear
    term beta_3*t represents long-term drift. H-IAR is subsequently fitted to
    the residuals rather than to the raw signal.
    """
    omega = 2 * np.pi / 365.25
    return beta_0 + beta_1 * np.cos(omega * t) + beta_2 * np.sin(omega * t) + beta_3 * t

def obter_parametros_sazonais(t, y):
    """
    Fit the deterministic seasonal/trend model to the training segment.

    Return the fitted parameters beta = (beta_0, beta_1, beta_2, beta_3).
    If nonlinear least-squares fitting raises RuntimeError, return None,
    allowing a later fallback to centring based only on the training mean.

    Reference: scipy.optimize.curve_fit.
    """
    # Initial values: beta_0 is the training mean;
    # annual harmonic and linear-trend coefficients start at zero.
    p0 = [np.mean(y), 0.0, 0.0, 0.0] 
    try:
        parametros_otimos, _ = curve_fit(modelo_fourier_tendencia, t, y, p0=p0)
        return parametros_otimos
    except RuntimeError:
        # Signal failure so the next step can apply the centring fallback.
        return None

# Common time reference: 2020-01-01 in days since the Unix epoch.
# The Phase B CSV stores 'time' in milliseconds; H-IAR uses days.
TEMPO_ZERO_DIAS = 1577836800000 / (1000 * 60 * 60 * 24)

def processar_pixel(args):
    """
    Process one pixel time series in a Phase B worker.

    Receive longitude, latitude and the pixel's table subset, then:
    1. reconstruct the irregular time axis;
    2. aggregate within-day observations;
    3. split chronologically into 90/10 training/test segments before deseasonalizing;
    4. fit the deterministic model using training data only;
    5. estimate Phi with H-IAR using training data only;
    6. evaluate the estimate through one-step-ahead prediction on the test segment.

    Splitting before deseasonalization prevents fitting on future test data:
    no deterministic or stochastic parameter is calibrated using test observations.
    """
    lon, lat, df_pixel = args
    
    try:
        # 1. Reconstruct time and merge within-day observations.
        df_pixel = df_pixel.copy()

        # Convert 'time' from milliseconds to days since 2020-01-01.
        df_pixel['tempo_dias'] = (df_pixel['time'] / (1000 * 60 * 60 * 24)) - TEMPO_ZERO_DIAS

        # Group observations from the same day to avoid orbital duplicates and near-zero Delta_t.
        df_pixel['dia_inteiro'] = np.floor(df_pixel['tempo_dias'])
        df_agregado = df_pixel.groupby('dia_inteiro').mean(numeric_only=True).reset_index()

        # Sort the aggregated observations before the temporal split.
        df_agregado = df_agregado.sort_values('tempo_dias')
        
        # 2. Split chronologically before fitting any deterministic component.
        # The first 90% calibrates the pipeline; the final 10% is withheld from fitting
        # for predictive evaluation. Temporal order is preserved.
        n_total = len(df_agregado)
        n_treino = int(n_total * 0.90)
        
        # Minimum information requirement for H-IAR calibration:
        # discard pixels with insufficient training observations after masking and aggregation.
        if n_treino < 30:
            return None

        # Chronological partition: training = observed past; test = withheld future.
        df_treino = df_agregado.iloc[:n_treino].copy()
        df_teste = df_agregado.iloc[n_treino:].copy()
        
        # 3. Remove deterministic components without fitting to test observations.
        # Fit each band's seasonality and trend using only the training data.
        # Apply the same fitted beta parameters to the test segment.
        for banda in ['B2', 'B3', 'B4', 'B8']:
            y_treino = df_treino[banda].values
            t_treino = df_treino['tempo_dias'].values
            
            p_opt = obter_parametros_sazonais(t_treino, y_treino)
            
            if p_opt is not None:
                # Remove the fitted deterministic component from the training data.
                df_treino[banda] = y_treino - modelo_fourier_tendencia(t_treino, *p_opt)
                # Apply the same coefficients to the test data without recalibration.
                t_teste = df_teste['tempo_dias'].values
                df_teste[banda] = df_teste[banda].values - modelo_fourier_tendencia(t_teste, *p_opt)
            else:
                # Leakage-free fallback: if the harmonic fit fails, subtract only the
                # training mean from both the training and test observations.
                media_treino = np.mean(y_treino)
                df_treino[banda] = y_treino - media_treino
                df_teste[banda] = df_teste[banda].values - media_treino

        # H-IAR input arrays:
        # X_* contains the four bands as quaternion-vector components;
        # t_* contains the corresponding irregular times.
        X_treino = df_treino[['B2', 'B3', 'B4', 'B8']].values
        t_treino = df_treino['tempo_dias'].values
        
        X_teste = df_teste[['B2', 'B3', 'B4', 'B8']].values
        t_teste = df_teste['tempo_dias'].values
        
        # 4. Estimate H-IAR parameters using training data only.
        # R_instrumental is on the percentage-reflectance scale of the exported bands:
        # 2.0 is an assumed standard deviation of 2 percentage points per component.
        R_instrumental = np.eye(4) * (2.0**2)

        # Calibrate Phi exclusively from X_treino and t_treino.
        res = estimar_hiar(X_treino, t_treino, R_in=R_instrumental)
        
        # Euclidean norm of the estimated quaternion, used as a temporal-persistence
        # descriptor; the output retains its historical resilience-related field name.
        phi_otimo = res['phi_otimo']
        norma_phi = np.linalg.norm(phi_otimo)
        
        # 5. One-step-ahead predictive evaluation on the test segment.
        # At each step F_{Delta_t} propagates the preceding observed residual;
        # errors are measured against the next observed residual, not a raw reflectance.
        erros_quadraticos = []
        rmse_b2, rmse_b3, rmse_b4, rmse_b8 = np.nan, np.nan, np.nan, np.nan
        
        # Evaluate forecasts only for Phi with norm below one.
        if norma_phi < 1.0: 
            F_tj_temp = np.zeros((4, 4), dtype=np.float64)
            
            # First test observation:
            # propagate the last observed training residual to the first test time.
            dt_inicial = t_teste[0] - t_treino[-1]
            calcular_matriz_transicao(phi_otimo, dt_inicial, F_tj_temp)
            x_pred = F_tj_temp @ X_treino[-1]
            erros_quadraticos.append((X_teste[0] - x_pred)**2)
            
            # Subsequent test observations:
            # condition each one-step forecast on the immediately preceding observed residual.
            # This evaluates local temporal transitions, not a
            # recursively simulated free-running forecast path.
            for i in range(1, len(t_teste)):
                dt = t_teste[i] - t_teste[i-1]
                calcular_matriz_transicao(phi_otimo, dt, F_tj_temp)
                
                # Predict observation i from the observed residual at i-1, not a filtered state.
                x_pred = F_tj_temp @ X_teste[i-1]
                erros_quadraticos.append((X_teste[i] - x_pred)**2)
                
            # Collect squared errors in an (N_teste, 4) array.
            # Each column corresponds to a band; compute RMSE separately
            # for each component.
            matriz_erros = np.array(erros_quadraticos)
            rmse_bandas = np.sqrt(np.mean(matriz_erros, axis=0))
            
            rmse_b2 = rmse_bandas[0]
            rmse_b3 = rmse_bandas[1]
            rmse_b4 = rmse_bandas[2]
            rmse_b8 = rmse_bandas[3]

        # Compatibility with result-key names from engine versions:
        # prefer 'iteracoes', accepting 'nit' if a SciPy-style key is supplied.
        iteracoes = res.get('iteracoes') if res.get('iteracoes') is not None else res.get('nit', 'N/A')
            
        # 6. Package the pixelwise spatial output.
        # The resulting table supports subsequent mapping:
        # Phi components, parameter norm, bandwise RMSE, numerical success and iterations.
        return {
            'Longitude': lon,
            'Latitude': lat,
            'N_Treino': len(t_treino),
            'N_Teste': len(t_teste),
            'Phi_a': phi_otimo[0],
            'Phi_b': phi_otimo[1],
            'Phi_c': phi_otimo[2],
            'Phi_d': phi_otimo[3],
            'Resiliencia_Global': norma_phi,
            'RMSE_B2': rmse_b2,
            'RMSE_B3': rmse_b3,
            'RMSE_B4': rmse_b4,
            'RMSE_B8': rmse_b8,
            'Convergencia': res['convergencia'],
            'Iteracoes': iteracoes
        }
    except Exception as e:
        # A local pixel exception does not stop regional processing.
        # Discard that pixel and continue parallel processing.
        return None

if __name__ == '__main__':
    # Main Phase B entry point: read the exported spatiotemporal table
    # and distribute independent pixelwise H-IAR tasks.
    print("Iniciando Fase B: Leitura da Malha Florestal Integral...")
    tempo_inicio = time.time()
    
    # Load the CSV produced by fase_b_dados.py:
    # longitude, latitude, B2, B3, B4, B8 and time.
    df_bruto = pd.read_csv('Extracao_MataSantaGenebra_Integral_Corrigida.csv')
    
    # Each unique coordinate pair defines a pixel time series.
    # Isolate pixels into tasks for independent parallel processing.
    tarefas = [(lon, lat, df_pixel) for (lon, lat), df_pixel in df_bruto.groupby(['longitude', 'latitude'])]
    total_pixeis = len(tarefas)
    
    print(f"Matriz Espacial Identificada: {total_pixeis} píxeis válidos.")
    print("Distribuindo processamento estocástico pelos núcleos do Servidor...")
    
    resultados_fase_B = []
    
    # Multicore process-based parallelism:
    # each worker receives an independent pixel series, avoiding shared series state
    # and accelerating regional estimation.
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futuros = {executor.submit(processar_pixel, tarefa): tarefa for tarefa in tarefas}
        
        # Collect results asynchronously as workers finish,
        # without imposing a spatial or temporal completion order.
        processados = 0
        for futuro in concurrent.futures.as_completed(futuros):
            processados += 1
            resultado = futuro.result()
            
            # Retain non-None worker results, including fits with convergencia=False.
            if resultado is not None:
                resultados_fase_B.append(resultado)

            # Report operational progress every 1000 pixels.
            if processados % 1000 == 0 or processados == total_pixeis:
                perc = (processados / total_pixeis) * 100
                print(f"Progresso: {processados}/{total_pixeis} ({perc:.1f}%) píxeis calculados.")
                
    print("Processamento concluído. Exportando termografia...")

    # Consolidate pixelwise results into a table for subsequent mapping.
    df_resultados_B = pd.DataFrame(resultados_fase_B)
    df_resultados_B.to_csv('Resultados_Termografia_SantaGenebra_FaseB.csv', index=False)
    
    # Compute total elapsed time in minutes.
    tempo_total = (time.time() - tempo_inicio) / 60
    print(f"Tarefa finalizada com sucesso! Tempo total: {tempo_total:.2f} minutos.")
