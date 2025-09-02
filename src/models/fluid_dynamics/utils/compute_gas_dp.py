import numpy as np
from math import sqrt

def compute_gas_pressure_drop(Qgas: float, p_inlet: float, dist: float, diam:float, density_GP:float, Tavg:float, T0: float, P0: float, pipe_eff: float, kw:float = 0.0037435):
    psi_to_atm = 0.0689476
    atm_to_mpa = 0.101325
    m3_to_cf = 35.3147

    # pf = np.sqrt((((qgas/1e3) ** 2)*dist) / ((ka**2) * (diam ** 5.33)) + pi**2)
    aux_expr =  ((1e3*Qgas/((m3_to_cf)*pipe_eff*(T0/(P0*1e3))*kw))/(1/sqrt(dist/1e3))*((diam*1e3)**2.667)**2)
    pf = sqrt(((p_inlet*psi_to_atm*atm_to_mpa*1e3)**2 - aux_expr*(density_GP*Tavg*((P0*1e3)/(0.375*T0))**2))/1e6)
    
    dp_gas = (p_inlet*psi_to_atm*atm_to_mpa - pf)
    
    return dp_gas