from .utils import *

def compute_multiphase_pressure_drop(pin: float):
    """
    parameters:
    pin: float | Presión en nodo inicial [MPa]
    dist: float | Distancia [???]
    """
    dp_gas = compute_gas_pressure_drop()
    dp_liq = compute_liquid_pressure_drop()
    ixlm = dp_gas/dp_liq
    ylp = ((ixlm**(1/4.12))+1)**4.12
    dp_mp = dp_liq*ylp
    return (dp_mp, ixlm, ylp)