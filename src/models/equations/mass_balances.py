from gamspy import Container, Domain, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

def mass_balances(m: Container) -> list[Equation]:
    n = m["n"]
    nn = m["nn"]
    d = m["d"]
    c = m["c"]
    t = m["t"]
    i = m["i"]
    j = m["j"]
    pf = m["pf"]
    arcs = m["arcs"]

    q_prod = m["Qprod"]
    q_inter = m["Qinter"]
    q_process = m["Qprocess"]
    
    mass_balance = Equation(m, "mass_balance", domain=[n,c,t], description="Mass balance per component at each node")
    mass_balance[i, c, t] =  q_prod[i, t, c] + Sum(Domain(nn, d).where[arcs[nn, i]], q_inter[nn, i, d, t, c]) == Sum(Domain(nn, d).where[arcs[i, nn]], q_inter[i, nn, d, t, c])
    mass_balance[j, c, t] =  Sum(Domain(nn, d).where[arcs[nn, j]], q_inter[nn, j, d, t, c]) == Sum(Domain(nn, d).where[arcs[j, nn]], q_inter[j, nn, d, t, c])
    mass_balance[pf, c, t] = Sum(Domain(nn, d).where[arcs[nn, pf]], q_inter[nn, pf, d, t, c]) == Sum(Domain(nn, d).where[arcs[pf, nn]], q_inter[pf, nn, d, t, c]) + q_process[pf, t, c]
    
    return [mass_balance]