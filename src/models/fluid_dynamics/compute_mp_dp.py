from .utils import *

def compute_multiphase_pressure_drop(Qoil, 
                                     Qgas, 
                                     Qwater: float, 
                                     p_inlet: float,
                                     dist: float,
                                     diam: float,
                                     density_LP: float = 955,
                                     viscosity_LP: float = 0.001,
                                     density_GP: float = 0.729,
                                     pipe_roughness: float = 0.0001,
                                     Tavg: float = 288.9,
                                     P0: float = 0.101325,
                                     T0: float = 288.9,
                                     n: float = 4.12,
                                     pipe_efficiency: float = 0.92):
    """
    parameters:
    Qoil: float | Flow production Oil [bbl/day]
    Qgas: float | Flow production Gas [kscf/day]
    Qwater: float | Flow production Water [bbl/day]
    p_inlet: float | Pressure at inlet node [PSI]
    distance: float | Distance [miles]
    diameter: float | Diameter [inches]
    density_LP: float | Density of the liquid phase [kg/m³]
    viscosity_LP: float | Viscosity of the liquid phase [Pa.s]
    density_GP: float | Density of the gas phase [kg/m³]
    pipe_roughness: float | Pipe internal wall roughness [m]
    Tavg: float | Average gas temperature [K]
    P0: float | Gas pressure at standard conditions [MPa]
    T0: float | Gas temperature at standard conditions [K]
    n: float | Flow regime of each phase, usually assumed turbulent
    pipe_efficiency: float | Pipe efficiency factor
    """
    inches_to_m = 0.0254
    mile_to_km = 1.60934
    mile_to_m = mile_to_km * 1000

    diam = diam * inches_to_m
    dist = dist * mile_to_m

    dp_gas = compute_gas_pressure_drop(Qgas, p_inlet, dist, diam, density_GP, Tavg, T0, P0, pipe_efficiency) # dp_gas [MPa/m]
    dp_liq = compute_liquid_pressure_drop(Qoil, Qwater, diam, density_LP, viscosity_LP, pipe_roughness)  # dp_liq [MPa]
    ixlm = dp_gas/dp_liq
    ylp = ((ixlm**(1/n))+1)**n
    dp_mp = dp_liq*ylp
    dp_mp = dp_mp*dist
    return (dp_mp, ixlm, ylp)