from gamspy import Container, Domain, Ord, Equation, Model, Sum, Sense, Options

def capacity_constraints(m: Container) -> list[Equation]:
    arcs = m["arcs"]
    pf = m["pf"]
    t = m["t"]
    tp = m["tp"]
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
    qGAS_interSQ = m["QGASinterSQ"]
    q_process = m["Qprocess"]
    accumulated_capacity = m["accumulated_capacity"]
    
    facility_capacity = Equation(m, "facility_capacity", domain=[pf,t,c], description="Facility must have enough processing capacity for each component")
    # facility_capacity[pf,t,c].where[tp[t]] = q_process[pf, t, c] <= Sum(Domain(s, tt).where[(Ord(tt) <= Ord(t)) & tp[tt]], (capacity[s,c])*y_pf[pf, s, tt]) # TODO: Daba problemas la misma expresión cuando tenía el Domain
    facility_capacity[pf,t,c].where[tp[t]] = q_process[pf, t, c] <= accumulated_capacity[pf,t,c]
    pipeline_capacity = Equation(m, "pipeline_capacity", domain=[n, nn, d, t, c], description="A pipeline must be installed to transport flow")
    pipeline_capacity[n,nn,d,t,c].where[arcs[n, nn] & tp[t]] = q_inter[n,nn,d,t,c] <= maxFlow[c,t]*Sum(tt.where[(Ord(tt) <= Ord(t)) & tp[tt]], x_bar[n,nn,d,tt])
    pipelineSQ_capacity = Equation(m, "pipelineSQ_capacity", domain=[n, nn, d, t], description="A pipeline must be installed to transport squared gas flow")
    pipelineSQ_capacity[n,nn,d,t].where[arcs[n, nn] & tp[t]] = qGAS_interSQ[n,nn,d,t] <= (maxFlow['gas',t]**2)*Sum(tt.where[(Ord(tt) <= Ord(t)) & tp[tt]], x_bar[n,nn,d,tt])
    # pipeline_capacity[n,nn,d,t,c].where[arcs[n, nn] & tp[t]] = q_inter[n,nn,d,t,c] <= maxFlow[c,t]*Sum(tt.where[Ord(tt) <= Ord(t)], x_bar[n,nn,d,tt])
    unique_pipeline = Equation(m, "unique_capacity", domain=[n], description="Just one pipeline can be installed for each connection")
    unique_pipeline[n] = Sum(Domain(nn,d,t).where[arcs[n, nn] & tp[t]], x_bar[n,nn,d,t]) <= 1
    define_acum_cap = Equation(m, "def_acum_cap", domain=[pf,c,t], description="Compute accumulated capacity of processing facility at each time period")
    define_acum_cap[pf,c,t] = Sum(Domain(s,tt).where[Ord(tt) <= Ord(t)], (capacity[s,c])*y_pf[pf, s, tt]) == accumulated_capacity[pf,t, c]

    # return [facility_capacity, pipeline_capacity, unique_pipeline]
    return [facility_capacity, pipeline_capacity, unique_pipeline, define_acum_cap]