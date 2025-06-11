from gamspy import Container, Domain, Ord, Equation, Model, Sum, Sense, Options

def capacity_constraints(m: Container) -> list[Equation]:
    pf = m["pf"]
    t = m["t"]
    tt = m["tt"]
    s = m["s"]
    n = m["n"]
    nn = m["nn"]
    d = m["d"]
    c = m["c"]
    
    capacity = m["capacity"]
    maxFlow = m["maxFlow"]

    y_pf = m["y_pf"]
    x_bar = m["x_bar"]
    q_inter = m["Qinter"]
    q_process = m["Qprocess"]
    
    facility_capacity = Equation(m, "facility_capacity", domain=[pf,t,c], description="Facility must have enough processing capacity for each component")
    facility_capacity[pf,t,c] = q_process[pf, t, c] <= Sum(Domain(s, tt).where[Ord(tt) <= Ord(t)], capacity[s,c]*y_pf[pf, s, tt])
    pipeline_capacity = Equation(m, "pipeline_capacity", domain=[n, nn, d, t, c], description="A pipeline must be installed to transport flow")
    pipeline_capacity[n,nn,d,t,c] = q_inter[n,nn,d,t,c] <= maxFlow[c,t]*Sum(tt.where[Ord(tt) <= Ord(t)], x_bar[n,nn,d,t])

    return [facility_capacity, pipeline_capacity]