"""
Monte Carlo evaluation of the H-IAR model
========================================
This script adapts the simulation-based evaluation approach used in the
IAR/CIAR/BIAR literature to the quaternion H-IAR model.

The experiment generates synthetic series at irregular observation times,
estimates Phi with the H-IAR engine, and aggregates absolute bias, standard
deviation, mean execution time and iteration counts by scenario and sample size.

Methodological references:
- Elorrieta et al. (2021): BIAR simulation design for irregular time series.
- Elorrieta, Eyheramendy & Palma (2019): state-space estimation, Kalman filtering
  and maximum likelihood.

Computational references:
- NumPy: pseudorandom generation and the exponential distribution.
- SciPy: L-BFGS-B optimization in the H-IAR estimator.
- Matplotlib: saving and explicitly closing figures.
"""

import numpy as np
import math
from numba import njit
from typing import Tuple
from scipy.optimize import minimize
from h_iar import calcular_matriz_transicao, estimar_hiar
import pandas as pd
from tqdm import tqdm
import time
import matplotlib.pyplot as plt

# ============================
# DATA-GENERATION FUNCTIONS
# ============================

def gerar_tempos_biar(N: int) -> np.ndarray:
    """
    Generate irregular observation times from a mixture of two exponential
    distributions, following the BIAR simulation design.

    Returns:
    - t: cumulative observation times, with t[0] = 0;
    - dt_array: irregular gaps Delta_t between consecutive observations.

    Reference: Elorrieta et al. (2021), BIAR Monte Carlo design.
    """
    # Scale parameters of the exponential mixture.
    # In np.random.exponential(scale=...), 'scale' is beta = 1/lambda.
    # The variable names lambda_1/lambda_2 follow the experiment's notation,
    # but their values are passed to NumPy as scales, not rates.
    lambda_1, lambda_2 = 15.0, 2.0
    w_1, w_2 = 0.15, 0.85
    
    dt_array = np.zeros(N - 1)
    for i in range(N - 1):
        # Select a mixture component: use the larger scale with probability w_1;
        # otherwise use the smaller scale. This reproduces heterogeneous
        # observation gaps in the BIAR simulation design.
        if np.random.rand() < w_1:
            dt_array[i] = np.random.exponential(scale=lambda_1)
        else:
            dt_array[i] = np.random.exponential(scale=lambda_2)
            
    # Set the first time to zero; subsequent times are cumulative sums
    # of Delta_t, giving an irregular observation grid.
    t = np.zeros(N)
    t[1:] = np.cumsum(dt_array)
    return t, dt_array

def gerar_serie_hiar(N: int, t: np.ndarray, dt_array: np.ndarray, phi_real: np.ndarray) -> np.ndarray:
    """
    Generate a synthetic H-IAR series under known parameters.

    The series is simulated in R^4, representing the four real components of
    a quaternion. At each irregular step, apply the matrix F_tj = Phi^Delta_tj
    and add Gaussian white noise with variance consistent with P0 = I_4.

    This routine is used only for Monte Carlo evaluation: the true parameter
    Phi is known and subsequently compared with its maximum likelihood estimate.
    """
    X = np.zeros((N, 4))
    
    # Isotropic initial condition: X_0 ~ N(0, I_4), consistent with P0 = I_4
    # in the Monte Carlo data-generating mechanism.
    X[0] = np.random.randn(4)
    
    norma_phi_sq = np.sum(phi_real**2)
    F_tj = np.zeros((4, 4))
    
    for i in range(1, N):
        dt = dt_array[i-1]
        
        # 1. Transition matrix for the current interval.
        calcular_matriz_transicao(phi_real, dt, F_tj)
        
        # 2. State-noise variance under P0 = I_4.
        # With isotropic generation, Q_tj = I_4*(1 - ||Phi||^(2*dt))
        # agrees with the fixed unconditional-covariance relation.
        variancia_ruido = 1.0 - (norma_phi_sq ** dt)
        
        # Isotropic Gaussian noise: epsilon_tj ~ N(0, variancia_ruido * I_4).
        # Generate four independent normal draws and multiply by the standard deviation.
        desvio_padrao_ruido = math.sqrt(max(variancia_ruido, 1e-12))
        epsilon_tj = np.random.randn(4) * desvio_padrao_ruido
        
        # 3. State transition without additional observation noise R in generation.
        X[i] = F_tj @ X[i-1] + epsilon_tj
        
    return X

# =====================================================================
# MONTE CARLO SIMULATION DRIVER
# =====================================================================

def executar_monte_carlo_hiar():
    """
    Run the Monte Carlo evaluation of the H-IAR estimator.

    The design follows the BIAR simulation approach: multiple parameter scenarios,
    different sample sizes and independent repetitions. For each combination,
    generate synthetic H-IAR series, estimate Phi, and summarize performance by:
    - mean parameter estimate;
    - absolute bias;
    - empirical standard deviation;
    - mean optimization time;
    - mean number of iterations.

    Reference: Elorrieta et al. (2021), Monte Carlo evaluation for irregular series.
    """
    # Fixed NumPy seed for the simulation draws; numerical results remain environment-dependent.
    np.random.seed(42)
    
    # Experimental grid: 1000 repetitions for each scenario and sample size.
    # N = 30, 100, 300 examines finite-sample behavior and
    # improvement with increasing sample size, not a formal asymptotic proof.
    repeticoes = 1000
    tamanhos_N = [30, 100, 300]
    
    # True parameter scenarios for Phi = (a,b,c,d).
    # The signs (+/+), (-/-), (-/+) and (+/-) combine scalar
    # and vector-component signs to examine estimation across
    # different quaternion orientations.
    casos_phi = {
        "Caso 1 (+/+)": np.array([0.7, 0.3, 0.3, 0.3]),
        "Caso 2 (-/-)": np.array([-0.7, -0.3, -0.3, -0.3]),
        "Caso 3 (-/+)": np.array([-0.9, 0.15, 0.15, 0.15]),
        "Caso 4 (+/-)": np.array([0.9, -0.15, -0.15, -0.15])
    }
    
    resultados_globais = []
    
    for nome_caso, phi_real in casos_phi.items():
        print(f"\nIniciando {nome_caso} | Phi Real: {phi_real}")
        
        for N in tamanhos_N:
            estimativas_phi = np.zeros((repeticoes, 4))
            tempos_execucao = np.zeros(repeticoes)
            iteracoes_array = np.zeros(repeticoes)
            
            loop_repeticoes = tqdm(range(repeticoes), desc=f"N={N}", leave=False)
            
            for m in loop_repeticoes:
                t, dt_array = gerar_tempos_biar(N)
                X_sim = gerar_serie_hiar(N, t, dt_array, phi_real)
                
                # Time the H-IAR estimation step only.
                # Exclude synthetic data generation from this timing to isolate
                # optimizer/filter cost.
                inicio = time.perf_counter()
                res = estimar_hiar(X_sim, t, R_in=np.eye(4)*1e-6)
                fim = time.perf_counter()
                
                estimativas_phi[m] = res['phi_otimo']
                tempos_execucao[m] = fim - inicio
                iteracoes_array[m] = res['iteracoes']
                
            # Parameterwise statistical summaries:
            # mean estimates, empirical bias, absolute bias and Monte Carlo
            # standard deviation (np.std uses ddof=0).
            media_estimada = np.mean(estimativas_phi, axis=0)
            bias = media_estimada - phi_real
            abs_bias = np.abs(bias)
            desvio_padrao = np.std(estimativas_phi, axis=0)
            
            # Computational summaries: mean estimation time, mean L-BFGS-B iterations
            # and mean time divided by the mean iteration count.
            tempo_medio = np.mean(tempos_execucao)
            iteracoes_medias = np.mean(iteracoes_array)
            # Guard the normalization if no iterations were recorded.
            tempo_por_iteracao_ms = (tempo_medio / max(iteracoes_medias, 1)) * 1000 
            
            for dim_idx, dim_name in enumerate(['a', 'b', 'c', 'd']):
                resultados_globais.append({
                    "Caso": nome_caso,
                    "N": N,
                    "Parâmetro": dim_name,
                    "Valor Real": phi_real[dim_idx],
                    "Média Estimada": media_estimada[dim_idx],
                    "Abs Bias": abs_bias[dim_idx],
                    "SD": desvio_padrao[dim_idx],
                    "Tempo Total (s)": tempo_medio,
                    "Iterações": iteracoes_medias,
                    "Tempo/Iteração (ms)": tempo_por_iteracao_ms
                })

    df_resultados = pd.DataFrame(resultados_globais)
    
    # Reshape only for terminal display:
    # full final metrics remain in df_resultados.
    print("\n\n=== TABELA DE RESULTADOS MONTE CARLO ===")
    tabela_pivot = df_resultados.pivot_table(
        index=["Caso", "N"], 
        columns=["Parâmetro"], 
        values=["Abs Bias", "SD", "Tempo Total (s)", "Iterações"]
    )
    print(tabela_pivot.round(4))
    
    return df_resultados

def gerar_graficos_ic(df_resultados: pd.DataFrame):
    """
    Generate diagnostic plots for the Monte Carlo experiment.

    The panels display, by sample size N:
    - absolute bias of the scalar component 'a';
    - empirical standard deviation of the scalar component 'a';
    - mean estimation time for each scenario.

    These plots visualize simulation summaries; they are not additional
    estimation steps.
    """
    # Plot scalar parameter 'a' as a representative to limit visual complexity.
    # Full a, b, c and d metrics remain in the CSV.
    df_plot = df_resultados[df_resultados['Parâmetro'] == 'a']
    
    fig, eixos = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("Validação Estocástica e Computacional - Modelo H-IAR", fontsize=16, weight='bold')

    # 1. Bias diagnostic: empirical |Bias| of scalar 'a' as N increases.
    for caso in df_plot['Caso'].unique():
        dados_caso = df_plot[df_plot['Caso'] == caso]
        eixos[0].plot(dados_caso['N'], dados_caso['Abs Bias'], marker='o', label=caso)
    eixos[0].set_title("Consistência Assintótica (Abs Bias vs N)")
    eixos[0].set_xlabel("Tamanho da Amostra (N)")
    eixos[0].set_ylabel("Módulo do Viés |Bias|")
    eixos[0].grid(True, linestyle='--', alpha=0.6)
    eixos[0].legend()

    # 2. Precision diagnostic: Monte Carlo dispersion of scalar 'a' as N increases.
    for caso in df_plot['Caso'].unique():
        dados_caso = df_plot[df_plot['Caso'] == caso]
        eixos[1].plot(dados_caso['N'], dados_caso['SD'], marker='s', label=caso)
    eixos[1].set_title("Ganhos de Precisão (SD vs N)")
    eixos[1].set_xlabel("Tamanho da Amostra (N)")
    eixos[1].set_ylabel("Desvio Padrão (SD)")
    eixos[1].grid(True, linestyle='--', alpha=0.6)
    eixos[1].legend()

    # 3. Computational diagnostic: mean estimation time by sample size.
    for caso in df_plot['Caso'].unique():
        dados_caso = df_plot[df_plot['Caso'] == caso]
        eixos[2].plot(dados_caso['N'], dados_caso['Tempo Total (s)'], marker='^', label=caso)
    eixos[2].set_title("Escalabilidade Computacional")
    eixos[2].set_xlabel("Tamanho da Amostra (N)")
    eixos[2].set_ylabel("Tempo Médio de Otimização (segundos)")
    eixos[2].grid(True, linestyle='--', alpha=0.6)
    eixos[2].legend()

    plt.tight_layout()
    plt.savefig("graficos_hiar_mc.png", dpi=300, bbox_inches='tight')
    plt.close()

# Run the full Monte Carlo experiment.
df_mc = executar_monte_carlo_hiar()

# Export the complete aggregate Monte Carlo table for reproducibility
# and subsequent external statistical analysis; individual fits are not saved.
df_mc.to_csv("resultados_hiar_mc.csv", index=False)

# Generate diagnostic plots.
gerar_graficos_ic(df_mc)
