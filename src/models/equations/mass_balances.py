from gamspy import Container, Domain, Ord, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

def mass_balances(m: Container) -> list[Equation]:
    n = m["n"]
    nn = m["nn"]
    d = m["d"]
    c = m["c"]
    t = m["t"]
    tt = m["tt"]
    i = m["i"]
    j = m["j"]
    pf = m["pf"]
    arcs = m["arcs"]

    q_prod = m["Qprod"]
    fluid_mult = m["fluid_mult"]
    st_time = m["st_time"]
    q_inter = m["Qinter"]
    q_process = m["Qprocess"]
    
    mass_balance_ij = Equation(m, "mass_balance_ij", domain=[i,c,t], description="Mass balance for source-junction")
    mass_balance_jpf = Equation(m, "mass_balance_jpf", domain=[j,c,t], description="Mass balance for junction-processing facility")
    mass_balance_pf = Equation(m, "mass_balance_pf", domain=[pf,c,t], description="Mass balance for processing facility")
    # mass_balance[i, c, t] =  q_prod[i, t, c] + Sum(Domain(nn, d).where[arcs[nn, i]], q_inter[nn, i, d, t, c]) == Sum(Domain(nn, d).where[arcs[i, nn]], q_inter[i, nn, d, t, c])
    # mass_balance[i, c, t] =  q_prod[i, t, c] == Sum(Domain(nn, d).where[arcs[i, nn]], q_inter[i, nn, d, t, c])
    mass_balance_ij[i, c, t].where[Ord(t) >= (st_time[i])] =  Sum(tt.where[Ord(tt) == (Ord(t) + 1 - st_time[i])], q_prod[c, tt]) == Sum(Domain(nn, d).where[arcs[i, nn]], q_inter[i, nn, d, t, c])
    mass_balance_jpf[j, c, t] =  Sum(Domain(nn, d).where[arcs[nn, j]], q_inter[nn, j, d, t, c]) == Sum(Domain(nn, d).where[arcs[j, nn]], q_inter[j, nn, d, t, c])
    # mass_balance[pf, c, t] = Sum(Domain(nn, d).where[arcs[nn, pf]], q_inter[nn, pf, d, t, c]) == Sum(Domain(nn, d).where[arcs[pf, nn]], q_inter[pf, nn, d, t, c]) + q_process[pf, t, c]
    mass_balance_pf[pf, c, t] = Sum(Domain(nn, d).where[arcs[nn, pf]], q_inter[nn, pf, d, t, c]) == q_process[pf, t, c]

    return [mass_balance_ij, mass_balance_jpf, mass_balance_pf]