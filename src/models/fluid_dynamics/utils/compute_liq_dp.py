from math import pi, log10, sqrt

def compute_liquid_pressure_drop(Qoil: float, Qwater: float, diam: float, density_LP: float, viscosity_LP:float, pipe_roughness:float):
    bbl_to_m3 = 0.158987
    day_to_sec = 86400

    velocity_liquid = ((Qoil + Qwater)/(pi*(diam/2)**2))*(bbl_to_m3/day_to_sec)
    reynolds_liquid = (density_LP * velocity_liquid * diam / viscosity_LP)/1E4
    
    aux_expr = pipe_roughness/(3.7*diam) + 6.9/(reynolds_liquid*1E4) # No se agrega el ^1.1 para ser más conservador
    friction_factor_liquid = (1/(-1.8*log10(aux_expr)))**2
    dp_liq = ((density_LP*(velocity_liquid**2)*friction_factor_liquid)/(2*diam))/1E6
    return dp_liq