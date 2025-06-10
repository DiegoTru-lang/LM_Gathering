from gamspy import Container, Domain, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

def costs_computation(m: Container) -> list[Equation]:
    n = m["n"]
    nn = m["nn"]
    d = m["d"]
    cpipe_km = m["cpipe_km"]
    t = m["t"]
    pipe_cost = m["pipe_cost"]
    x_bar = m["x_bar"]
    dist = m["dist"]
    arcs = m["arcs"]

    compute_pipe_cost = Equation(m, "compute_pipe_cost_per_t", domain=t, description="Compute total pipeline installation cost at period 't' ")
    compute_pipe_cost[t] =  pipe_cost[t] == Sum(d, cpipe_km[d]*Sum(Domain(n,nn).where[arcs[n, nn]], dist[n,nn]*x_bar[n,nn,d,t]))
    
    return [compute_pipe_cost]