"""
H-IAR core engine
=================
This module implements parameter estimation for H-IAR, a quaternion
extension of the IAR/CIAR/BIAR family for irregularly observed time series.

Each observation is represented by a four-dimensional real vector
associated with the quaternion q_t = a_t + b_t i + c_t j + d_t k.

The irregular transition uses Phi^Delta_t, evaluated through the polar
representation and real powers of quaternions, and converted into a real
4x4 matrix of left multiplication. The near-real branch follows the
convention documented in calcular_matriz_transicao and docs/METHOD_NOTES.md.

Main methodological references:
- Elorrieta, Eyheramendy & Palma (2019): CIAR, state-space representation,
  weak stationarity and maximum likelihood through the Kalman filter.
- Elorrieta, Eyheramendy, Palma & Ojeda (2021): BIAR and multicomponent
  autoregressive structure in irregular time series.
- Kuipers (1999): polar representation and real powers in H*.
- Mebius (2005): matrix representations of quaternion multiplication in R^4.
"""

import numpy as np
import math
from numba import njit
from typing import Tuple
from scipy.optimize import minimize

@njit
def calcular_matriz_transicao(phi_params: np.ndarray, dt: float, F_tj: np.ndarray) -> np.ndarray:
    """
    Build the 4x4 transition matrix associated with a temporal quaternion power.

    Parameters:
    phi_params (np.ndarray): Components [a, b, c, d] of quaternion Phi.
    dt (float): Interval between consecutive observations.
    F_tj (np.ndarray): Preallocated 4x4 array, filled in place.

    Returns:
    np.ndarray: The same F_tj array, representing left multiplication by Phi^dt.

    For a vector-part norm below 1e-12, the implemented scalar power is a**dt
    when a > 0, and zero otherwise; all vector coefficients are zero.
    """
    a, b, c, d = phi_params
    
    # Vector part Im(Phi) = (b, c, d), defining the unit axis in polar form.
    norma_v_sq = b**2 + c**2 + d**2
    norma_v = math.sqrt(norma_v_sq)
    
    # Degenerate polar case: Im(Phi) near zero prevents a unique definition of axis u.
    # For positive real Phi, the power reduces to a^dt; nonpositive a uses zero here.
    if norma_v < 1e-12:
        A = (a)**dt if a > 0 else 0.0 
        B, C, D = 0.0, 0.0, 0.0
    else:
        # Compute the total norm using the already computed vector part.
        norma_phi = math.sqrt(a**2 + norma_v_sq)

        # Polar form of Phi: Phi = ||Phi|| (cos(theta) + u sin(theta)).
        # Reference: Kuipers (1999), polar representation and real powers in H*.
        theta = math.acos(a / norma_phi)
        
        # Unit axis u = Im(Phi)/||Im(Phi)|| in the polar representation.
        u_b, u_c, u_d = b / norma_v, c / norma_v, d / norma_v
        
        # Real power Phi^dt:
        # radial scaling ||Phi||^dt and angular advance theta -> dt*theta.
        escala_tempo = norma_phi**dt
        angulo_tempo = theta * dt
        
        A = escala_tempo * math.cos(angulo_tempo)
        seno_tempo = escala_tempo * math.sin(angulo_tempo)
        
        B = u_b * seno_tempo
        C = u_c * seno_tempo
        D = u_d * seno_tempo

    # Real 4x4 matrix for left multiplication by Phi^dt = A + Bi + Cj + Dk.
    # This operator represents the Hamilton product on the real vector (a,b,c,d).
    # Reference: Mebius (2005), quaternion multiplication matrices in R^4.
    F_tj[0, 0] = A;  F_tj[0, 1] = -B; F_tj[0, 2] = -C; F_tj[0, 3] = -D
    F_tj[1, 0] = B;  F_tj[1, 1] = A;  F_tj[1, 2] = -D; F_tj[1, 3] = C
    F_tj[2, 0] = C;  F_tj[2, 1] = D;  F_tj[2, 2] = A;  F_tj[2, 3] = -B
    F_tj[3, 0] = D;  F_tj[3, 1] = -C; F_tj[3, 2] = B;  F_tj[3, 3] = A
    
    return F_tj

@njit
def filtro_kalman_hiar(
    X: np.ndarray, 
    delta_t: np.ndarray, 
    phi_params: np.ndarray, 
    norma_phi_sq: float,
    R: np.ndarray, 
    X0: np.ndarray, 
    P0: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Run the Kalman-filter recursion for the H-IAR state-space model.

    Parameters:
    X (np.ndarray): Observation matrix with shape (N, 4).
    delta_t (np.ndarray): Precomputed time intervals with shape (N-1,).
    phi_params (np.ndarray): Components [a, b, c, d] of quaternion Phi.
    norma_phi_sq (float): Squared norm of quaternion Phi.
    R (np.ndarray): Observation-error covariance matrix with shape (4, 4).
    X0 (np.ndarray): Initial state vector with shape (4,).
    P0 (np.ndarray): Initial state covariance matrix with shape (4, 4).

    Returns:
    v (np.ndarray): Innovation vectors with shape (N, 4).
    Lambda_inv_mat (np.ndarray): Inverse innovation covariances, shape (N, 4, 4).
    det_Lambda_vec (np.ndarray): Innovation covariance determinants, shape (N,).

    Index zero is reserved for the initial condition; recursion starts at j=1.
    """
    N = X.shape[0]
    
    # Preallocate arrays reused throughout the recursion.
    v = np.zeros((N, 4), dtype=np.float64)
    Lambda_inv_mat = np.zeros((N, 4, 4), dtype=np.float64)
    det_Lambda_vec = np.zeros(N, dtype=np.float64)
    F_tj = np.zeros((4, 4), dtype=np.float64)
    
    # Identity matrix shared by regularization and update steps.
    I4 = np.eye(4, dtype=np.float64)
    
    # Initial state and covariance.
    X_hat = X0.copy()
    P = P0.copy()
    
    # Recursion begins at j=1; j=0 is reserved for the anchored initial condition.
    for j in range(1, N):
        dt = delta_t[j-1]
        
        # 1. Dynamic transition matrix.
        calcular_matriz_transicao(phi_params, dt, F_tj)
        
        # 2. State-noise covariance under the approximate diagonal specification.
        # The exact fixed-covariance relation would be Q_j = P0 - F_j P0 F_j.T.
        # Here Q_j = P0 * (1 - ||Phi||^(2*dt)) uses the empirical marginal
        # variances as scale factors and is diagonal nonnegative for admissible norms.
        # This does not generally preserve P0 as a stationary covariance when anisotropic.
        #
        # This is the adopted H-IAR specification for anisotropic data:
        # it does not impose equal spectral variances and avoids indefinite Q_j
        # from the fixed-covariance expression under quaternion rotations.
        Q_tj = P0 * (1.0 - norma_phi_sq**dt)
        
        # Prediction step.
        X_pred = F_tj @ X_hat
        P_pred = F_tj @ P @ F_tj.T + Q_tj
        
        # Observation and innovation step.
        # All four components are observed, so G = I_4.
        # Thus v = X - G*X_pred simplifies to v = X - X_pred.
        v_tj = X[j] - X_pred
        
        # Innovation covariance: Lambda = G*P_pred*G.T + R -> P_pred + R.
        Lambda_tj = P_pred + R
        
        det_L = np.linalg.det(Lambda_tj)
        
        # Numerical stabilization through a diagonal perturbation (jitter) of Lambda_tj.
        # When det(Lambda_tj) reaches the numerical threshold, add 1e-6 * I_4
        # before evaluating log(det) and the matrix inverse.
        if det_L <= 1e-12:
            Lambda_tj = Lambda_tj + (I4 * 1e-6)
            det_L = np.linalg.det(Lambda_tj)
        
        # Store quantities reused in the likelihood calculation.
        v[j] = v_tj
        det_Lambda_vec[j] = det_L
        
        # Update step.
        # Kalman gain: K = P_pred * G.T * Lambda^-1 -> P_pred * Lambda^-1.
        # Store the inverse because it is already needed for the Kalman gain
        # and will be reused in the NLL calculation.
        Lambda_inv = np.linalg.inv(Lambda_tj)
        Lambda_inv_mat[j] = Lambda_inv
        
        K_tj = P_pred @ Lambda_inv
        
        # Update the state and covariance for the next step.
        X_hat = X_pred + K_tj @ v_tj
        P = (I4 - K_tj) @ P_pred

        # Correct loss of symmetry caused by floating-point rounding.
        P = (P + P.T) * 0.5
        
    return v, Lambda_inv_mat, det_Lambda_vec

@njit
def calcular_nll_hiar(v: np.ndarray, Lambda_inv_mat: np.ndarray, det_Lambda_vec: np.ndarray) -> float:
    """
    Compute the Gaussian negative log-likelihood from Kalman-filter innovations.

    For each j >= 1, accumulate:
    log|Lambda_j| + v_j.T @ Lambda_j^{-1} @ v_j + k log(2*pi),
    where k = 4 observed components. Return one half of the accumulated sum.

    Reference: the innovation likelihood formulation for state-space models,
    as used in CIAR for maximum likelihood estimation through the Kalman filter.
    """
    N_total = v.shape[0]
    
    # Observation dimension: four real quaternion components.
    k = 4.0 
    
    # Exclude the initial slot from the sum because it represents the anchor.
    nll_soma = (N_total - 1) * (k * np.log(2.0 * np.pi))
    
    for j in range(1, N_total):
        v_tj = v[j]
        
        # Retrieve the determinant and inverse already computed by the filter.
        det_Lambda = det_Lambda_vec[j]
        Lambda_inv = Lambda_inv_mat[j]
            
        # Gaussian quadratic term: v_j.T @ Lambda_j^{-1} @ v_j.
        termo_quadratico = v_tj @ Lambda_inv @ v_tj
        
        # Accumulate the terms depending on Lambda_j and v_j.
        nll_soma += np.log(det_Lambda) + termo_quadratico
        
    return 0.5 * nll_soma

@njit
def funcao_objetivo(phi_params_flat: np.ndarray, X: np.ndarray, delta_t: np.ndarray, 
                    R: np.ndarray, X0: np.ndarray, P0: np.ndarray) -> float:
    """
    Evaluate the objective minimized by L-BFGS-B: the Gaussian Kalman-innovation
    NLL with radial projection and a quadratic penalty outside the radial limit.
    """
    # Ensure contiguous memory layout for SciPy/Numba interoperability.
    phi_params_flat = np.ascontiguousarray(phi_params_flat)
    
    norma_sq = phi_params_flat[0]**2 + phi_params_flat[1]**2 + \
               phi_params_flat[2]**2 + phi_params_flat[3]**2

    limite_estacionariedade = 0.99
    
    # 1. Radial projection and continuous quadratic penalty outside the norm limit.
    # L-BFGS-B bounds define a componentwise box; likelihood evaluation here uses
    # projection into ||Phi||^2 <= 0.99 plus a penalty for the unprojected excess.
    if norma_sq > limite_estacionariedade:
        fator_escala = math.sqrt(limite_estacionariedade / norma_sq)
        phi_valido = phi_params_flat * fator_escala
        norma_sq_valida = limite_estacionariedade
        
        # Anchor the excess exactly at the boundary to avoid a jump in the penalty.
        excesso = norma_sq - limite_estacionariedade
        penalidade = 1e6 * (excesso ** 2)
    else:
        phi_valido = phi_params_flat
        norma_sq_valida = norma_sq
        penalidade = 0.0
    
    # 2. Kalman recursion with admissible parameters.
    v, Lambda_inv_mat, det_Lambda_vec = filtro_kalman_hiar(
        X, delta_t, phi_valido, norma_sq_valida, R, X0, P0
    )
    
    # 3. Evaluate the negative log-likelihood.
    nll = calcular_nll_hiar(v, Lambda_inv_mat, det_Lambda_vec)
    
    return nll + penalidade

def estimar_hiar(X: np.ndarray, t: np.ndarray, R_in: np.ndarray = None) -> dict:
    """
    Estimate the parameters of the H-IAR model.
    """
    # Preprocessing: componentwise centring.
    # Remove the empirical mean to estimate H-IAR dynamics around X0 = 0.
    X = X - np.mean(X, axis=0)
    
    # --- STEP 1: PREPROCESSING AND INITIAL VALUES ---
    
    phi_init = np.array([0.5, 0.0, 0.0, 0.0], dtype=np.float64)
    
    # Empirical diagonal P0: use each component's empirical variance.
    # This retains unequal spectral scales and anchors the approximate
    # diagonal specification used for Q_tj.
    variancias_marginais = np.var(X, axis=0)
    P0 = np.diag(variancias_marginais)
    
    # Observation/instrument-error covariance matrix.
    # Use a small observation covariance if R_in is not supplied.
    if R_in is None:
        R = np.eye(4, dtype=np.float64) * 1e-6
    else:
        R = R_in
        
    # Zero initial state.
    X0 = np.zeros(4, dtype=np.float64)
    
    # Precompute irregular intervals Delta_t = t_j - t_{j-1}.
    delta_t = np.diff(t)
    
    # --- STEP 2: NUMERICAL OPTIMIZATION ---
    
    # Componentwise bounds; likelihood evaluation uses the radial limit ||Phi||^2 <= 0.99
    # through projection plus a penalty inside the objective function.
    limites = [(-0.99, 0.99), (-0.99, 0.99), (-0.99, 0.99), (-0.99, 0.99)]
    
    resultado_otimizacao = minimize(
        fun=funcao_objetivo,
        x0=phi_init,
        args=(X, delta_t, R, X0, P0),
        method='L-BFGS-B',
        bounds=limites,
        options={'disp': True, 'maxiter': 2000, 'ftol': 1e-9} 
    )
    
    phi_otimo = resultado_otimizacao.x
    
    # --- STEP 3: POST-FIT QUANTITIES ---
    
    # Final radial projection into ||Phi||^2 <= 0.99
    # if the soft penalty allowed termination outside that limit.
    norma_otima_sq = np.sum(phi_otimo**2)
    if norma_otima_sq > 0.99:
        phi_otimo = phi_otimo * np.sqrt(0.99 / norma_otima_sq)
        norma_otima_sq = 0.99
    
    # Innovations associated with the estimated parameter.
    v_final, _, _ = filtro_kalman_hiar(
        X, delta_t, phi_otimo, norma_otima_sq, R, X0, P0
    )
    
    # Discard the slot reserved for the initial condition.
    v_validos = v_final[1:]
    
    # Base state-noise matrix for Delta_t = 1,
    # consistent with the approximate diagonal specification used for Q_tj.
    Sigma_otimo = P0 * (1.0 - norma_otima_sq)
    
    # Post-fit diagnostic of contemporaneous correlation among filtered innovations.
    # Rho_xi does not enter the NLL; it summarizes residual intercomponent dependence.
    Rho_xi = np.corrcoef(v_validos, rowvar=False)
    
    resultados = {
        'phi_otimo': phi_otimo,
        'nll_minima': resultado_otimizacao.fun,
        'Sigma_otimo': Sigma_otimo, 
        'Rho_xi': Rho_xi,
        'convergencia': resultado_otimizacao.success,
        'mensagem': resultado_otimizacao.message,
        'iteracoes': resultado_otimizacao.nit,
        'nfev': resultado_otimizacao.nfev
    }
    
    return resultados
