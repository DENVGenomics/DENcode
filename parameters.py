
import random
import math
import numpy as np

# === Parameter Sampler ===
def sample_params(seed=None):
    if seed is not None:
        random.seed(seed)
    return {
        "T_i": random.uniform(25, 30),
        "beta_h": random.uniform(0.7, 0.8),
        "P_mh": random.uniform(0.7, 0.8),
        "delta": random.uniform(0.006, 0.009),
    }

# === Temperature-dependent biological functions ===

def viremia_curve(t):
    # Normalized Gaussian decay: peak at 2 days, rapid decline
    return np.exp(-0.5 * ((t - 2) / 1.5) ** 2)

def eip_temperature(T, beta_0=2.9, beta_T=-0.08, tau=4.9):
    return math.exp(beta_0 + beta_T * T + 0.5 * tau)

def mosquito_mortality(T):
    return 0.8692 - 0.159 * T + 0.01116 * T**2 - 3.408e-4 * T**3 + 3.809e-6 * T**4

def mosquito_survival(T):
    mu = mosquito_mortality(T)
    eip = eip_temperature(T)
    return math.exp(-mu * eip)
