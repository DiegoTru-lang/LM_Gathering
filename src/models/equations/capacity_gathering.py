from gamspy import Container, Domain, Ord, Equation, Model, Sum, Sense, Options

def capacity_constraints(m: Container) -> list[Equation]:
    arcs = m["arcs"]
    pf = m["pf"]
    t = m["t"]
    tp = m["tp"]
    tt = m["tt"]
    n = m["n"]
    nn = m["nn"]
    d = m["d"]
    c = m["c"]
    
    maxFlow = m["maxFlow"]

    x_bar = m["x_bar"]
    q_inter = m["Qinter"]
    
    pipeline_capacity = Equation(m, "pipeline_capacity", domain=[n, nn, d, t, c], description="A pipeline must be installed to transport flow")
    pipeline_capacity[n,nn,d,t,c].where[arcs[n, nn] & tp[t]] = q_inter[n,nn,d,t,c] <= maxFlow[c]*Sum(tt.where[(Ord(tt) <= Ord(t)) & tp[tt]], x_bar[n,nn,d,tt])

    
    unique_pipeline = Equation(m, "unique_capacity", domain=[n], description="Just one pipeline can be installed for each connection")
    unique_pipeline[n] = Sum(Domain(nn,d,t).where[arcs[n, nn] & tp[t]], x_bar[n,nn,d,t]) <= 1

    return [pipeline_capacity, unique_pipeline]