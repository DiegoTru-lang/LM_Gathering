from math import sqrt

def compute_weymouth_constant(dist: float, diam: float, density_GP: float = 0.729, Tavg: float = 288.9, T0: float = 288.9, P0: float = 0.101325, pipe_eff: float = 0.92, kw: float = 0.0037435, z: float = 1.0) -> float:
    """
    Units:
    distance: float | Distance [miles]
    diameter: float | Diameter [inches]
    density_GP: float | Density of the gas phase [kg/m³]
    Tavg: float | Average gas temperature [K]
    P0: float | Gas pressure at standard conditions [MPa]
    T0: float | Gas temperature at standard conditions [K]
    pipe_efficiency: float | Pipe efficiency factor
    """
    inches_to_m = 0.0254
    mile_to_km = 1.60934
    
    dist = dist * mile_to_km # miles to km
    diam = diam * inches_to_m * 1e3 # inches to mm

    wey_const =  (sqrt(density_GP*Tavg*(dist)*z))/(kw*pipe_eff*(T0/(P0*1e3))*((diam)**2.667))

    return wey_const