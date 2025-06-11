from gamspy import Container, Domain, Ord, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

def costs_computation(m: Container) -> list[Equation]:
    n = m["n"]
    nn = m["nn"]
    pf = m["pf"]
    s = m["s"]
    d = m["d"]
    t = m["t"]
    dist = m["dist"]
    arcs = m["arcs"]

    cpipe_km = m["cpipe_km"]
    facility_cost_size = m["facility_cost_size"]
    ir = m["ir"]

    x_bar = m["x_bar"]
    y_pf = m["y_pf"]
    pipe_cost = m["pipe_cost"]
    facility_cost = m["facility_cost"]

    total_cost = m["total_cost"]

    compute_pipe_cost = Equation(m, "compute_pipe_cost_per_t", domain=t, description="Compute total pipeline installation cost at period 't' ")
    compute_pipe_cost[t] =  pipe_cost[t] == Sum(d, cpipe_km[d]*Sum(Domain(n,nn).where[arcs[n, nn]], dist[n,nn]*x_bar[n,nn,d,t]))
    
    compute_facility_cost = Equation(m, "compute_facility_cost_per_t", domain=t, description="Compute total facility installation cost at period 't' ")
    compute_facility_cost[t] = facility_cost[t] == Sum([pf,s], facility_cost_size[s]*y_pf[pf, s, t])

    compute_total_cost = Equation(m, "compute_total_cost", description="Compute total discounted cost")
    compute_total_cost[...] = total_cost == Sum(t, (1 + ir)**(-(Ord(t) - 1)) * (pipe_cost[t] + facility_cost[t]))

    return [compute_pipe_cost, compute_facility_cost, compute_total_cost]