import numpy as np
from math import sqrt

def compute_gas_pressure_drop(Qgas: float, p_inlet: float, dist: float, diam:float, density_GP:float, Tavg:float, T0: float, P0: float, pipe_eff: float, kw:float = 0.0037435, z:float = 1.0):
    psi_to_atm = 0.0689476
    psi_to_kpa = 6.89476
    psi_to_mpa = 0.00689476
    atm_to_mpa = 0.101325
    m3_to_cf = 35.3147

    # pf = np.sqrt((((qgas/1e3) ** 2)*dist) / ((ka**2) * (diam ** 5.33)) + pi**2)
    # aux_expr =  (((1e3*Qgas/(m3_to_cf))/(pipe_eff*(T0/(P0*1e3))*kw))/((1/sqrt(dist/1e3))*((diam*1e3)**2.667)))**2
    # pf = sqrt(((p_inlet*psi_to_atm*atm_to_mpa*1e3)**2 - aux_expr*(density_GP*Tavg*((P0*1e3)/(0.375*T0))**2))/1e6)
    aux_expr =  ((Qgas*1e3/m3_to_cf)*sqrt(density_GP*Tavg*(dist/1e3)*z))/(kw*pipe_eff*(T0/(P0*1e3))*((diam*1e3)**2.667))
    pf = sqrt(max((p_inlet*psi_to_kpa)**2 - aux_expr**2, 0))
    pf = pf/1e3 #kPa to MPa
    
    dp_gas = (p_inlet*psi_to_mpa - pf)/dist
    # dp_gas = (p_inlet*psi_to_mpa - pf)

    return dp_gas