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
    x_bar = m["x_bar"]
    selected_pipes = m["sel_pipes"]
    dp_max = m["dp_max"]
    deltaP = m["deltaP"]
    deltaPgas = m["deltaPgas"]
    deltaPliq = m["deltaPliq"]
    deltaPliqYLM = m["deltaPliqYLM"]
    ixlm_ub = m["ixlm_ub"]
    ylp = m["ylp"]
    allowed_int = m["allowed_int"]

    one_interval = Equation(m, "one_interval", domain=[j,pf,d,t], description= "Each connection (j, pf, d) must lie between a single 'pw' at each time period 't' ")
    one_interval[j,pf,d,t].where[selected_pipes[j,pf,d] & tp[t]] =  Sum(pw.where[Ord(pw) <= allowed_int[j,pf,d,t]], x_pw[j,pf,d,t,pw]) == Sum(tt.where[Ord(tt) <= Ord(t)], x_bar[j,pf,d,tt])
    # one_interval[j,pf,d,t].where[selected_pipes[j,pf,d] & tp[t]] =  Sum(pw.where[Ord(pw) <= allowed_int[j,pf,d]], x_pw[j,pf,d,t,pw]) == 1
    
    ixlm_interval = Equation(m, "ixlm_interval", domain=[j,pf,d,t,pw], description= "Impose interval for LM correlation between junction 'j' and processing facility 'pf' during time period 't' ")
    ixlm_interval[j,pf,d,t,pw].where[selected_pipes[j,pf,d] & tp[t] & (Ord(pw) <= allowed_int[j,pf,d,t])] =  deltaPgas[j,pf,t] <= (ixlm_ub[pw,j,pf,d,t] - 1e3)*deltaPliq[j,pf,d,t] + dp_max*(1 - x_pw[j,pf,d,t,pw])
    # ixlm_interval[j,pf,d,t,pw].where[selected_pipes[j,pf,d] & tp[t] & (Ord(pw) <= allowed_int[j,pf,d])] =  deltaPgas[j,pf,t] <= ixlm_ub[pw,j,pf,d]*Sum(d, deltaPliq[j,pf,d,t]) + dp_max*(1 - x_pw[j,pf,t,pw])

    compute_deltaPliqYLM = Equation(m, "compute_deltaPliqYLM", domain=[j,pf,d,t], description="Intermediate variable definition")
    compute_deltaPliqYLM[j,pf,d,t].where[selected_pipes[j,pf,d] & tp[t]] =  deltaPliq[j,pf,d,t] == Sum(pw.where[Ord(pw) <= allowed_int[j,pf,d,t]], deltaPliqYLM[j,pf,d,t,pw])
    # compute_deltaPliqYLM[j,pf,t].where[selected_pipes[j,pf,d] & tp[t]] =  Sum(d, deltaPliq[j,pf,d,t]) == Sum(pw.where[Ord(pw) <= allowed_int[j,pf,d]], deltaPliqYLM[j,pf,t,pw])
    
    deltaPliqYLM_bound = Equation(m, "deltaPliqYLM_bound", domain=[j,pf,d,t,pw], description="Fix to 0 intermediate variable if IXLM does not lie between interval 'pw' ")
    deltaPliqYLM_bound[j,pf,d,t,pw].where[selected_pipes[j,pf,d] & tp[t] & (Ord(pw) <= allowed_int[j,pf,d,t])] = deltaPliqYLM[j,pf,d,t,pw] <= dp_max*x_pw[j,pf,d,t,pw]

    pressure_drop_from_j_to_pf = Equation(m, "pressure_drop_from_j_to_pf", domain=[j,pf,d,t], description= "Compute pressure drop from junction 'j' to processing facility 'pf' ")
    pressure_drop_from_j_to_pf[j,pf,d,t].where[selected_pipes[j,pf,d] & tp[t]] =  deltaP[j,pf,t] >= Sum(pw.where[Ord(pw) <= allowed_int[j,pf,d,t]], ylp[pw,j,pf,d,t]*deltaPliqYLM[j,pf,d,t,pw])

    return [one_interval, ixlm_interval, compute_deltaPliqYLM, deltaPliqYLM_bound, pressure_drop_from_j_to_pf]