from gamspy import Container, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

def mass_balances(m: Container):
    """
    supply = Equation(m, "supply", domain=m.i, description="...",)
    supply[m.i] = Sum(m.j, m.x[m.i, m.j] <= m.cap)
    demand = Equation(m, "demand", domain=m.j, description="...",)
    demand[m.j] = Sum(m.i, m.x[m.i, m.j] >= m.dem)
    """
    pass