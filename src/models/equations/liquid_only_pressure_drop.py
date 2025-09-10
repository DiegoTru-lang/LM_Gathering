from gamspy import Container, Domain, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options
from gamspy.math import sqrt
from numpy import pi

def liquid_only_pressure_drop(m: Container) -> list[Equation]:
    j = m["j"]
    pf = m["pf"]
    d = m["d"]
    t = m["t"]
    tp = m["tp"]
    arcs = m["arcs"]
    deltaPliq = m["deltaPliq"]
    hffl = m["hffl"]
    rho_liq = m["rho_liq"]
    vel_liq = m["vel_liq"]
    dist = m["dist"]
    diam = m["diam"]
    q_inter = m["Qinter"]

    in_to_m = 0.0254
    mile_to_km = 1.60934
    bbl_to_m3 = 0.158987
    compute_velLIQ = Equation(m, "compute_velLIQ", domain=[j,pf,d,t], description= "Compute liquid velocity assuming liquid-only on pipeline ")
    compute_velLIQ[j,pf,d,t].where[arcs[j,pf]  & tp[t]] =  vel_liq[j, pf, d, t] >= ((q_inter[j, pf, d, t, "oil"] + q_inter[j, pf, d, t, "water"])*bbl_to_m3)/((pi*((diam[d]*in_to_m)**2)/4)*24*3600)

    compute_dpLIQ = Equation(m, "compute_dpLIQ", domain=[j,pf,d,t], description= "Compute liquid-only pressure drop between nodes 'j' and 'pf' during time period 't' [MPa]")
    compute_dpLIQ[j,pf,d,t].where[arcs[j,pf] & tp[t]] =  deltaPliq[j, pf, d, t]*1e3 >= (dist[j,pf]*mile_to_km)*(hffl/2)*(rho_liq*vel_liq[j,pf,d,t]*vel_liq[j,pf,d,t])/(diam[d]*in_to_m) # RHS in kPa

    return [compute_velLIQ, compute_dpLIQ]