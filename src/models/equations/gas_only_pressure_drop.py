from gamspy import Container, Domain, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options
from gamspy.math import sqrt

def gas_only_pressure_drop(m: Container) -> list[Equation]:
    j = m["j"]
    pf = m["pf"]
    d = m["d"]
    t = m["t"]
    arcs = m["arcs"]
    selected_pipes = m["sel_pipes"]
    deltaPgas = m["deltaPgas"]
    press = m["press"]
    pressSQ = m["pressSQ"]
    pressGAS = m["pressGAS"]
    qGAS_interSQ = m["QGASinterSQ"]
    kw = m["kw"]

    m3_to_cf = 35.3147

    weymouth_correlation = Equation(m, "weymouth_correlation", domain=[j,pf,d,t], description= "Weymouth correlation for gas-only pressure drop between nodes 'n' and 'nn' during time period 't' assuming diameter 'd' [MPa]")
    # weymouth_correlation = Equation(m, "weymouth_correlation", domain=[n,nn,t], description=" Weymouth correlation for gas-only pressure drop between nodes 'n' and 'nn' during time period 't' [MPa]")
    # Units: pressSQ [MPa^2] --> [kPa^2]; qGAS_interSQ [(kscf/d)^2] -> [(m3/d)^2]
    # weymouth_correlation[j,pf,d,t].where[arcs[j,pf]] =  pressGAS[pf, t] <= sqrt(pressSQ[j,t]*(1e3**2) - (qGAS_interSQ[j, pf, d, t]*((1e3/m3_to_cf)**2))*kw[j,pf,d])/1e3 # kPa --> MPa
    weymouth_correlation[j,pf,d,t].where[selected_pipes[j,pf,d]] =  (pressGAS[pf, t]*pressGAS[pf, t])*(1e3**2) <= pressSQ[j,t]*(1e3**2) - (qGAS_interSQ[j, pf, d, t]*((1e3/m3_to_cf)**2))*kw[j,pf,d] # kPa --> MPa

    compute_dpGAS = Equation(m, "compute_dpGAS", domain=[j,pf,t], description= "Compute gas-only pressure drop between nodes 'j' and 'pf' during time period 't' [MPa]")
    compute_dpGAS[j,pf,t] =  deltaPgas[j, pf, t] == press[j,t] - pressGAS[pf,t]
    # compute_dpGAS[j,pf,t].where[selected_pipes[j,pf,d]] =  deltaPgas[j, pf, t] == press[j,t] - pressGAS[pf,t]

    return [weymouth_correlation, compute_dpGAS]