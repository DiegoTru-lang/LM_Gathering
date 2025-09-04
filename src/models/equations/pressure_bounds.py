from gamspy import Container, Domain, Set, Parameter, Ord, Variable, Equation, Model, Sum, Sense, Options

def pressure_bounds(m: Container) -> list[Equation]:
    i = m["i"]
    j = m["j"]
    pf = m["pf"]
    d = m["d"]
    t = m["t"]
    tp = m["tp"]
    tt = m["tt"]
    x_bar = m["x_bar"]
    arcs = m["arcs"]
    press = m["press"]
    pressGAS = m["pressGAS"]
    pressSQ = m["pressSQ"]
    deltaP = m["deltaP"]
    deltaPgas = m["deltaPgas"]
    fixPress = m["fixPress"]
    maxPress = m["maxPress"]
    pmin_pf = m["pmin_pf"]

    compute_pressure_junction = Equation(m, "compute_pressure_junction", domain=[i,j,d,t], description= "Compute max pressure at junction 'j' based on connections from sources ")
    compute_pressure_junction[i,j,d,t].where[arcs[i,j] & tp[t]] =  press[j,t] <= fixPress[i,j,d,t] + maxPress*(1 - Sum(tt.where[Ord(tt) <= Ord(t)], x_bar[i,j,d,tt]))

    compute_square_pressure_junction = Equation(m, "compute_square_pressure_junction", domain=[i,j,d,t], description= "Compute max squared pressure at junction 'j' based on connections from sources ")
    compute_square_pressure_junction[i,j,d,t].where[arcs[i,j] & tp[t]] =  pressSQ[j,t] <= fixPress[i,j,d,t]*fixPress[i,j,d,t] + maxPress*maxPress*(1 - Sum(tt.where[Ord(tt) <= Ord(t)], x_bar[i,j,d,tt]))

    minimum_pressure_pf = Equation(m, "min_press_pf", domain=[pf,t], description= "Impose minimum pressure at processing facility 'pf' during time period 't'")
    minimum_pressure_pf[pf,t] =  press[pf,t] >= pmin_pf

    minimum_pressure_pf_GAS = Equation(m, "min_press_pf_GAS", domain=[pf,t], description= "Impose minimum pressure at processing facility 'pf' during time period 't'")
    minimum_pressure_pf_GAS[pf,t] =  pressGAS[pf,t] >= pmin_pf

    pressure_at_pf = Equation(m, "press_at_pf", domain=[j,pf,d,t], description= "Compute pressure at processing facility 'pf' during time period 't'")
    pressure_at_pf[j,pf,d,t].where[arcs[j,pf]] =  press[pf,t] <= press[j,t] - deltaP[j,pf,d,t] + maxPress*(1 - Sum(tt.where[Ord(tt) <= Ord(t)], x_bar[j,pf,d,tt]))

    return [compute_pressure_junction, compute_square_pressure_junction, minimum_pressure_pf, minimum_pressure_pf_GAS, pressure_at_pf]
    # return [minimum_pressure_pf, minimum_pressure_pf_GAS]