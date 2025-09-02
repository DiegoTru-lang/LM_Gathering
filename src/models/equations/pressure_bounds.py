from gamspy import Container, Domain, Set, Parameter, Ord, Variable, Equation, Model, Sum, Sense, Options

def gas_only_pressure_drop(m: Container) -> list[Equation]:
    i = m["i"]
    j = m["j"]
    pf = m["pf"]
    d = m["d"]
    t = m["t"]
    tt = m["tt"]
    x_bar = m["x_bar"]
    arcs = m["arcs"]
    press = m["press"]
    pressGAS = m["pressGAS"]
    pressSQ = m["pressSQ"]
    fixPress = m["fixPress"]
    maxPress = m["maxPress"]
    pmin_pf = m["pmin_pf"]

    compute_pressure_junction = Equation(m, "compute_pressure_junction", domain=[i,j,t], description= "Compute max pressure at junction 'j' based on connections from sources ")
    compute_pressure_junction[i,j,t].where[arcs[i,j]] =  press[j,t] <= fixPress[i,j,t] + maxPress*(1 - Sum(Domain(d,tt).where[Ord(tt) <= Ord(t)], x_bar[i,j,d,t]))

    compute_square_pressure_junction = Equation(m, "compute_square_pressure_junction", domain=[i,j,t], description= "Compute max squared pressure at junction 'j' based on connections from sources ")
    compute_square_pressure_junction[i,j,t].where[arcs[i,j]] =  pressSQ[j,t] <= fixPress[i,j,t]*fixPress[i,j,t] + maxPress*maxPress*(1 - Sum(Domain(d,tt).where[Ord(tt) <= Ord(t)], x_bar[i,j,d,t]))

    minimum_pressure_pf = Equation(m, "min_press_pf", domain=[pf,t], description= "Impose minimum pressure at processing facility 'pf' during time period 't'")
    minimum_pressure_pf[pf,t] =  press[pf,t] >= pmin_pf

    minimum_pressure_pf_GAS = Equation(m, "min_press_pf_GAS", domain=[pf,t], description= "Impose minimum pressure at processing facility 'pf' during time period 't'")
    minimum_pressure_pf_GAS[pf,t] =  pressGAS[pf,t] >= pmin_pf

    return [compute_pressure_junction, compute_square_pressure_junction, minimum_pressure_pf, minimum_pressure_pf_GAS]