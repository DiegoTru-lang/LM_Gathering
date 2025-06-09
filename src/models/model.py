from .data_classes import DataClass
from gamspy import Container, Set, Parameter, Variable, Equation, Model, Sum, Sense, Options

class Model():
    def __init__(self, model_name: str, data: DataClass = None):
        self.model_name = model_name
        self.m = _construct_model()
        if data:
            self.instance_model(data)

    def _construct_model(self, data: DataClass) -> Container:
        m = Container()

        n = Set(m, "n", description="Nodes in the network")
        nn = Set(m, "nn", domain=n)
        i = Set(m, "i", domain=n, description="Source nodes / multi-phase")
        j = Set(m, "j", domain=n, description="Junction nodes / multi-phase")
        pf = Set(m, "pf", domain=n, description="Processing facilities")

        arcs = Set(m, "arcs", domain=[n,nn], description="Allowed connections") 
        
        t = Set(m, "t", description="Time periods")
        c = Set(m, "c", description="Flow components")
        d = Set(m, "d", description="Pipeline diameter options")
        s = Set(m, "s", description="Facility sizes")

        prod = Parameter(m, "prod", domain=[i,t,c], description="Production at source node 'i' of component 'c' during time period 't' [mscf per day]")
        st_time = Parameter(m, "st_time", domain=i, description="Production start time of source node 'i' ")
        capacity = Parameter(m, "capacity", domain=[s, c], description="Capacity for facility size 's' and component 'c' [mscf per day]")
        facility_cost = Parameter(m, "facility_cost", domain=s, description="Cost for facility size 's' [$]")
        # install_facility_cost = Parameter(m, "install_facility_cost", domain=pf, description="Cost for installing a facility in node 'pf' [$]")
        cpipe = Parameter(m, "cpipe", description="Pipeline cost per inch per km [$/(in*km)]")
        cpipe_km = Parameter(m, "cpipe_km", domain=d, description="Pipeline cost per km for diameter size 'd' [$/km]")
        mult_pipe = Parameter(m, "mult_pipe", domain=d, description="Multiplier for pipeline of diameter 'd'")
        diam = Parameter(m, "diam", domain=d, description="Diameter of pipeline option 'd' [inches]")
        ir = Parameter(m, "ir", description="Interest rate [%]")

        loc_x = Parameter(m, "loc_x", domain=n, description="Position of node 'n' in x-axis [km]")
        loc_y = Parameter(m, "loc_y", domain=n, description="Position of node 'n' in y-axis [km]")
        dist = Parameter(m, "dist", domain=[n,nn], description="Distance between nodes'n'to 'nn' [km]")

        pwell = Parameter(m, "pwell", domain=[i,t], description="Wellhead pressure per source node 'i' at time period 't' [MPa]")
        pmin_pf = Parameter(m, "pmin_pf", description="Minimum inlet pressure at processing facility [MPa]")
        
        # These variables must instantiate just for [n, nn] included in arcs set 
        y_pf = Variable(m, "y_pf", domain=[pf, s, t], type="binary", 
                        description="Equals 1 if a processing facility of size 's' is installed at node 'pf' during time period 't'")
        x = Variable(m, "x", domain=[n, nn, d, t], type="binary", 
                     description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        q = Variable(m, "Q", domain=[n, nn, d, t, c], description="Flow of component 'c' through pipeline segment between nodes 'n' and 'nn' of diameter 'd' during time period 't' [mscf per day]")
        
        # add topology constraints topology_constraints(m)
        # list_eqs = topology_constraints(m)
        # add rest of constraints

        # add objective function

        # model = Model(m, "Two-echelon Multiphase Pipeline network design model",
        #               equations=[...]", # eqs = list_eqs + list_eqs2
        #               sense=Sense.MIN, options=Options.SOLVER="cplex")
        #               objective=obj,
        #               options=Options.SOLVER="cplex")
        return m

    def instance_model(self, data: DataClass):
        self.m.i.set_records()
    
    def solve(self):
        ...

    def predict(self, input_data):
        print(f"Predicting with {self.model_name} for input: {input_data}")
        return "prediction_result"