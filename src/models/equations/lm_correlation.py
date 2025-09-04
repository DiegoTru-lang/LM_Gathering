from gamspy import Container, Domain, Set, Parameter, Ord, Variable, Equation, Model, Sum, Sense, Options

def lm_correlation(m: Container) -> list[Equation]:
    j = m["j"]
    pf = m["pf"]
    d = m["d"]
    t = m["t"]
    tp = m["tp"]
    tt = m["tt"]
    pw = m["pw"]
    x_pw = m["x_pw"]
    arcs = m["arcs"]
    dp_max = m["dp_max"]
    deltaP = m["deltaP"]
    deltaPgas = m["deltaPgas"]
    deltaPliq = m["deltaPliq"]
    deltaPliqYLM = m["deltaPliqYLM"]
    ixlm_ub = m["ixlm_ub"]
    ylp = m["ylp"]

    one_interval = Equation(m, "one_interval", domain=[j,pf,t], description= "Each connection (j, pf) must lie between a single 'pw' at each time period 't' ")
    one_interval[j,pf,t].where[arcs[j,pf] & tp[t]] =  Sum(pw, x_pw[j,pf,t, pw]) == 1
    
    ixlm_interval = Equation(m, "ixlm_interval", domain=[j,pf,t, pw], description= "Impose interval for LM correlation between junction 'j' and processing facility 'pf' during time period 't' ")
    ixlm_interval[j,pf,t, pw].where[arcs[j,pf] & tp[t]] =  deltaPgas[j,pf,t] <= ixlm_ub[pw]*Sum(d, deltaPliq[j,pf,d,t]) + dp_max*(1 - x_pw[j,pf,t, pw])

    compute_deltaPliqYLM = Equation(m, "compute_deltaPliqYLM", domain=[j,pf,t], description="Intermediate variable definition")
    compute_deltaPliqYLM[j,pf,t].where[arcs[j,pf] & tp[t]] =  Sum(d, deltaPliq[j,pf,d,t]) == Sum(pw, deltaPliqYLM[j,pf,t,pw])
    
    deltaPliqYLM_bound = Equation(m, "deltaPliqYLM_bound", domain=[j,pf,t,pw], description="Fix to 0 intermediate variable if IXLM does not lie between interval 'pw' ")
    deltaPliqYLM_bound[j,pf,t,pw].where[arcs[j,pf] & tp[t]] = deltaPliqYLM[j,pf,t,pw] <= dp_max*x_pw[j,pf,t, pw]

    pressure_drop_from_j_to_pf = Equation(m, "pressure_drop_from_j_to_pf", domain=[j,pf,t], description= "Compute pressure drop from junction 'j' to processing facility 'pf' ")
    pressure_drop_from_j_to_pf[j,pf,t].where[arcs[j,pf] & tp[t]] =  deltaP[j,pf,t] >= Sum(pw, ylp[pw]*deltaPliqYLM[j,pf,t,pw])

    return [one_interval, ixlm_interval, compute_deltaPliqYLM, deltaPliqYLM_bound, pressure_drop_from_j_to_pf]