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

    x_conn = m["x_conn"]
    q_prod = m["Qprod"]
    st_time = m["st_time"]
    q_inter = m["Qinter"]
    q_processed = m["Qprocessed"]

    total_processed = m["total_processed"]
    
    mass_balance_ij = Equation(m, "mass_balance_ij", domain=[i,j,c,t], description="Mass balance for source-junction")
    mass_balance_jpf = Equation(m, "mass_balance_jpf", domain=[j,pf, c,t], description="Mass balance for junction-processing facility")
    compute_q_processed = Equation(m, "mass_balance_pf", domain=[c,t], description="Mass balance for junction-processing facility")


    compute_total_processed = Equation(m, "compute_total_processed", description="Compute total production")

    mass_balance_ij[i, j, c, t].where[Ord(t) >= (st_time[i])] =  Sum(tt.where[Ord(tt) == (Ord(t) + 1 - st_time[i])], q_prod[i, tt, c])*x_conn[i] == Sum(d, q_inter[i, j, d, t, c])
    
    mass_balance_jpf[j, pf, c, t] =  Sum(Domain(i, d), q_inter[i, j, d, t, c]) == Sum(d, q_inter[j, pf, d, t, c])
    compute_q_processed[c, t] =  Sum(i.where[Ord(t) >= (st_time[i])], x_conn[i]*Sum(tt.where[Ord(tt) == (Ord(t) + 1 - st_time[i])], q_prod[i, tt, c])) == q_processed[t, c]

    compute_total_processed[...] = total_processed == Sum(Domain(t, c), q_processed[t, c])

    return [mass_balance_ij, mass_balance_jpf, compute_q_processed, compute_total_processed]