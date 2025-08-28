import numpy as np

def compute_gas_pressure_drop(qgas: float, pi: float, dist: float, diam:float, ka:float = 5430.86):
    """
    Unidades:
    qgas -> Gas flow []
    pi -> Init pressure [MPa]
    dist -> Distance [km]
    diam -> Diameter [m]
    ka -> Constant for Weymouth
    """
    pf = np.sqrt((((qgas/1e3) ** 2)*dist) / ((ka**2) * (diam ** 5.33)) + pi**2)
    return pi - pf