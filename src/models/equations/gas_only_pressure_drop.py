from gamspy import Container, Domain, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

def gas_only_pressure_drop(m: Container) -> list[Equation]:
    n = m["n"]
    nn = m["nn"]
    d = m["d"]
    t = m["t"]
    x_bar = m["x_bar"]
    dist = m["dist"]
    arcs = m["arcs"]
    deltaPgas = m["deltaPgas"]
    press = m["press"]
    cpipe_km = m["cpipe_km"]

    weymouth_correlation = Equation(m, "weymouth_correlation", domain=[n,nn,t], description=" Weymouth correlation for gas-only pressure drop between nodes 'n' and 'nn' during time period 't' [MPa]")
    weymouth_correlation[n,nn,t].where[arcs[n,nn]] =  deltaPgas[n, nn, t] == press[n,t] - Sum(d, cpipe_km[d]*Sum(Domain(n,nn).where[arcs[n, nn]], dist[n,nn]*x_bar[n,nn,d,t]))
    
    return [...]