from gamspy import Container, Domain, Set, Parameter, Ord, Variable, Equation, Model, Sum, Sense, Options

def lm_correlation(m: Container) -> list[Equation]:
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
    deltaP = m["deltaP"]
    deltaPgas = m["deltaPgas"]
    deltaPliq = m["deltaPliq"]

    pressure_drop_from_j_to_pf = Equation(m, "pressure_drop_from_j_to_pf", domain=[j,pf,d,t], description= "Compute pressure drop from junction 'j' to processing facility 'pf' ")
    pressure_drop_from_j_to_pf[j,pf,d,t].where[arcs[j,pf] & tp[t]] =  deltaP[j,pf,d,t] >= deltaPliq[j,pf,d,t]

    return [pressure_drop_from_j_to_pf]