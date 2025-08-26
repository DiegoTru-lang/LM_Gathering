from .data_classes import DataClass
from .data_classes.utils.model_mappings import ModelSetMapping
from gamspy import Container, Set, Sum, Parameter, Variable, Alias, Model, Sense, Options, Problem
from .equations import mass_balances, costs_computation, capacity_constraints

class GatheringModel():
    def __init__(self, model_name: str, data: DataClass = None):
        self.model_name = model_name
        self.m = self._construct_model()
        if data:
            self.instance_model(data)

    def _construct_model(self) -> Container:
        m = Container()

        # Sets

        n = Set(m, "n", description="Nodes in the network")
        nn = Alias(m, name="nn", alias_with=n)
        i = Set(m, "i", domain=n, description="Source nodes / multi-phase")
        j = Set(m, "j", domain=n, description="Junction nodes / multi-phase")
        pf = Set(m, "pf", domain=n, description="Processing facilities")

        arcs = Set(m, "arcs", domain=[n,nn], description="Allowed connections") 
        
        t = Set(m, "t", description="Time periods")
        tt = Alias(m, name="tt", alias_with=t)
        c = Set(m, "c", records=["oil", "gas", "water"], description="Flow components")
        d = Set(m, "d", description="Pipeline diameter options")
        s = Set(m, "s", description="Facility sizes")

        # Parameters

        q_prod = Parameter(m, "Qprod", domain=t, description="Oil production at each node at 't' periods after start-time [BBL per day]")
        # q_prod = Parameter(m, "Qprod", domain=[i,t,c], description="Production at source node 'i' of component 'c' during time period 't' [mscf per day]")
        st_time = Parameter(m, "st_time", domain=i, description="Production start time of source node 'i' ") #TODO: Compute based on q_prod
        capacity = Parameter(m, "capacity", domain=[s, c], description="Capacity for facility size 's' and component 'c' [mscf per day]")
        diam = Parameter(m, "diam", domain=d, description="Diameter of pipeline option 'd' [inches]")
        facility_cost_size = Parameter(m, "facility_cost_size", domain=s, description="Cost for facility size 's' [$]")
        diameter_cost = Parameter(m, "diameter_cost", domain=d, description="Pipeline cost per km for diameter size 'd' [$/mile]")
        # install_facility_cost = Parameter(m, "install_facility_cost", domain=pf, description="Cost for installing a facility in node 'pf' [$]")
        # cpipe = Parameter(m, "cpipe", description="Pipeline cost per inch per km [$/(in*km)]")
        # mult_pipe = Parameter(m, "mult_pipe", domain=d, description="Multiplier for pipeline of diameter 'd'")
        # cpipe_km = Parameter(m, "cpipe_km", domain=d, description="Pipeline cost per km for diameter size 'd' [$/km]")
        # cpipe_km[d] = cpipe * mult_pipe[d] * diam[d]

        ir = Parameter(m, "ir", description="Interest rate [%]")

        loc_x = Parameter(m, "loc_x", domain=n, description="Position of node 'n' in x-axis [km]")
        loc_y = Parameter(m, "loc_y", domain=n, description="Position of node 'n' in y-axis [km]")
        dist = Parameter(m, "dist", domain=[n,nn], description="Distance between nodes'n' to 'nn' [mile]")

        fluid_mult = Parameter(m, "fluid_mult", domain=c, description="Multiplier to convert oil production to each component")
        fluid_mult["oil"] = 1
        fluid_mult["gas"] = 2
        fluid_mult["water"] = 3.5
        maxFlow = Parameter(m, "maxFlow", domain=[c,t], description="Maximum flow of component 'c' during time period 't' [mscf per day]") # TODO: Compute as sum of production at all source nodes
        maxFlow[c,t] = Sum(i, q_prod[t])*fluid_mult[c]
        # pwell = Parameter(m, "pwell", domain=[i,t], description="Wellhead pressure per source node 'i' at time period 't' [MPa]")
        pwell = Parameter(m, "pwell", domain=t, description="Wellhead pressure at every node at time period 't' [MPa]")
        pmin_pf = Parameter(m, "pmin_pf", description="Minimum inlet pressure at processing facility [MPa]")
        
        # Variables

        # These variables must instantiate just for [n, nn] included in arcs set 
        y_pf = Variable(m, "y_pf", domain=[pf, s, t], type="binary", 
                        description="Equals 1 if a processing facility of size 's' is installed at node 'pf' during time period 't'")
        x_bar = Variable(m, "x_bar", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        # xf = Variable(m, "xf", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' transports fluid from node 'n' to 'nn' at time period 't'")
        # xr = Variable(m, "xr", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' transports fluid from node 'nn' to 'n' at time period 't'")
        # xnf = Variable(m, "xnf", domain=[n, nn, d, t], type="binary", description="Equals 1 if no pipeline segment of diameter 'd' is installed between 'n' and 'nn' and P(n) <= P(nn) at time period 't'")
        # xnr = Variable(m, "xnr", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        
        q_inter = Variable(m, "Qinter", domain=[n, nn, d, t, c], description="Flow of component 'c' through pipeline segment between nodes 'n' and 'nn' of diameter 'd' during time period 't' [mscf per day]")
        q_process = Variable(m, "Qprocess", domain=[pf, t, c], description="Amount of component 'c' processed at facility 'pf' during time period 't' [mscf per day]")
        press = Variable(m, "press", domain=[n, t], description="Pressure at node 'n' during time period 't' [MPa]")
        deltaP = Variable(m, "deltaP", domain=[n, nn, t], description="Pressure drop 'multiphase' between nodes 'n' and 'nn' during time period 't' [MPa]")
        deltaPgas = Variable(m, "deltaPgas", domain=[n, nn, t], description="Pressure drop 'gas-only' between nodes 'n' and 'nn' during time period 't' [MPa]")
        deltaPliq = Variable(m, "deltaPliq", domain=[n, nn, t], description="Pressure drop 'liquid-only' between nodes 'n' and 'nn' during time period 't' [MPa]")

        pipe_cost = Variable(m, "pipe_cost", domain=t, description="Total cost on pipeline installation during time period 't' [$]")
        facility_cost = Variable(m, "facility_cost", domain=t, description="Total cost on facility installation during time period 't' [$]")
        total_cost = Variable(m, "total_cost", description="Total discounted cost [$]")
        
        # Call equations
        eqs = []
        eqs += mass_balances(m)
        eqs += costs_computation(m)
        eqs += capacity_constraints(m)

        model = Model(
            m,
            "Multiphase_network_design",
            # equations=eqs,
            equations=m.getEquations(),
            sense=Sense.MIN,
            problem=Problem.MIP,
            objective=total_cost,
            # options=Options(mip="gurobi"),
        )

        return m

    def instance_model(self, data: dict):
        # Automatically set records for all sets/parameters defined in ModelSetMapping
        for mapping in ModelSetMapping:
            model_attr, data_key = mapping.value
            if type(data_key) is list:
                combined_data = []
                for key in data_key:
                    if key in data:
                        if isinstance(data[key], list):
                            combined_data.extend(data[key])
                        else:
                            combined_data.append(data[key])
                self.m[model_attr].setRecords(combined_data)
                continue
            self.m[model_attr].setRecords(data.get(data_key))


    def compute_first_echelon_parameters(self):
        ...

    def solve(self):
        ...

    def predict(self, input_data):
        print(f"Predicting with {self.model_name} for input: {input_data}")
        return "prediction_result"