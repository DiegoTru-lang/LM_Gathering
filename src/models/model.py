from .data_classes import DataClass
from .data_classes.utils.model_mappings import ModelSetMapping
from gamspy import Container, Set, Sum, Parameter, Variable, Alias, Model, Sense, Ord, Options, Problem
from .equations import *
from .fluid_dynamics import compute_multiphase_pressure_drop, compute_weymouth_constant
from itertools import product
import sys

class GatheringModel():
    def __init__(self, model_name: str, data: DataClass = None):
        self.model_name = model_name
        self.m = self._construct_model()
        if data:
            self.data = data
            self.instance_model(data)

    def _construct_model(self) -> Container:
        m = Container(output=sys.stdout)

        # Sets

        n = Set(m, "n", description="Nodes in the network")
        nn = Alias(m, name="nn", alias_with=n)
        i = Set(m, "i", domain=n, description="Source nodes / multi-phase")
        j = Set(m, "j", domain=n, description="Junction nodes / multi-phase")
        pf = Set(m, "pf", domain=n, description="Processing facilities")

        arcs = Set(m, "arcs", domain=[n,nn], description="Allowed connections") 
        
        t = Set(m, "t", description="Time periods")
        tt = Alias(m, name="tt", alias_with=t)
        tp = Set(m, "tp", domain=t,description="Time periods where there is a production peak")
        c = Set(m, "c", records=["oil", "gas", "water"], description="Flow components")
        d = Set(m, "d", description="Pipeline diameter options")
        s = Set(m, "s", description="Facility sizes")

        # Parameters

        q_prod = Parameter(m, "Qprod", domain=[c,t], description="Oil production at each node at 't' periods after start-time [BBL per day]")
        # q_prod = Parameter(m, "Qprod", domain=[i,t,c], description="Production at source node 'i' of component 'c' during time period 't' [mscf per day]")
        st_time = Parameter(m, "st_time", domain=i, description="Production start time of source node 'i' ") #TODO: Compute based on q_prod
        capacity = Parameter(m, "capacity", domain=[s, c], description="Capacity for facility size 's' and component 'c' [mscf per day]")
        diam = Parameter(m, "diam", domain=d, description="Diameter of pipeline option 'd' [inches]")
        facility_cost_size = Parameter(m, "facility_cost_size", domain=s, description="Cost for facility size 's' [MMUSD]")
        diameter_cost = Parameter(m, "diameter_cost", domain=d, description="Pipeline cost per km for diameter size 'd' [kUSD/mile]")
        # install_facility_cost = Parameter(m, "install_facility_cost", domain=pf, description="Cost for installing a facility in node 'pf' [$]")
        # cpipe = Parameter(m, "cpipe", description="Pipeline cost per inch per km [$/(in*km)]")
        # mult_pipe = Parameter(m, "mult_pipe", domain=d, description="Multiplier for pipeline of diameter 'd'")
        # cpipe_km = Parameter(m, "cpipe_km", domain=d, description="Pipeline cost per km for diameter size 'd' [$/km]")
        # cpipe_km[d] = cpipe * mult_pipe[d] * diam[d]

        ir = Parameter(m, "ir", description="Interest rate [%]")

        loc_x = Parameter(m, "loc_x", domain=n, description="Position of node 'n' in x-axis [km]")
        loc_y = Parameter(m, "loc_y", domain=n, description="Position of node 'n' in y-axis [km]")
        dist = Parameter(m, "dist", domain=[n,nn], description="Distance between nodes'n' to 'nn' [mile]")

        hffl = Parameter(m, "hffl", records=0.025, description="Hydraulic friction factor for liquid phase")
        rho_liq = Parameter(m, "rho_liq", records=850, description="Density of liquid phase mixture")
        kw = Parameter(m, "kw", domain=[n,nn,d], description="Weymouth constants and parameters synthesized")
        maxFlow = Parameter(m, "maxFlow", domain=[c,t], description="Maximum flow of component 'c' during time period 't' [mscf per day]") # TODO: Compute as sum of production at all source nodes
        # pwell = Parameter(m, "pwell", domain=[i,t], description="Wellhead pressure per source node 'i' at time period 't' [MPa]")
        dp_firstEch = Parameter(m, "dp_firstEch", domain=[i,j,d,t], description="Pre-computed pressure drop at first echelon for connection (i,j), diameter 'd' at time period 't' [MPa]")
        dp_max = Parameter(m, "dp_max", description="Maximum allowable pressure drop between two nodes [MPa]")
        pwell = Parameter(m, "pwell", domain=t, description="Wellhead pressure at every node at time period 't' [MPa]")
        pmin_pf = Parameter(m, "pmin_pf", description="Minimum inlet pressure at processing facility [MPa]")
        # fixPress = Parameter(m, "fixPress", domain=[n, nn, t], description="Pre-computed max pressure at junction 'j' during time period 't' if connection (i, j) is installed [MPa]")
        fixPress = Parameter(m, "fixPress", domain=[i, j, d, t], description="Pre-computed max pressure at junction 'j' during time period 't' if connection (i, j) is installed with diameter 'd' [MPa]")
        maxPress = Parameter(m, "maxPress", description="Max pressure at node [MPa]")
        
        # Variables

        # These variables must instantiate just for [n, nn] included in arcs set 
        y_pf = Variable(m, "y_pf", domain=[pf, s, t], type="binary", 
                        description="Equals 1 if a processing facility of size 's' is installed at node 'pf' during time period 't'")
        x_bar = Variable(m, "x_bar", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        # xf = Variable(m, "xf", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' transports fluid from node 'n' to 'nn' at time period 't'")
        # xr = Variable(m, "xr", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' transports fluid from node 'nn' to 'n' at time period 't'")
        # xnf = Variable(m, "xnf", domain=[n, nn, d, t], type="binary", description="Equals 1 if no pipeline segment of diameter 'd' is installed between 'n' and 'nn' and P(n) <= P(nn) at time period 't'")
        # xnr = Variable(m, "xnr", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        
        q_inter = Variable(m, "Qinter", type="positive", domain=[n, nn, d, t, c], description="Flow of component 'c' through pipeline segment between nodes 'n' and 'nn' of diameter 'd' during time period 't' [mscf per day]")
        qGAS_interSQ = Variable(m, "QGASinterSQ", type="positive", domain=[n, nn, d, t], description="Squared flow of 'gas' through pipeline segment between nodes 'n' and 'nn' of diameter 'd' during time period 't' [(mscf per day)**2]")
        q_process = Variable(m, "Qprocess", type="positive", domain=[pf, t, c], description="Amount of component 'c' processed at facility 'pf' during time period 't' [mscf per day]")
        press = Variable(m, "press", type="positive", domain=[n, t], description="Pressure at node 'n' during time period 't' [MPa]")
        pressSQ = Variable(m, "pressSQ", type="positive", domain=[n, t], description="Squared pressure at node 'n' during time period 't' [MPa]")
        # pressGAS = Variable(m, "pressGAS", type="positive", domain=[n, t], description="Pressure at node 'n' during time period 't' assuming gas-only pressure drop [MPa]")
        pressGAS = Variable(m, "pressGAS", type="positive", domain=[pf, t], description="Pressure at node 'pf' during time period 't' assuming gas-only pressure drop [MPa]")
        deltaP = Variable(m, "deltaP", type="positive", domain=[n, nn, d, t], description="Pressure drop 'multiphase' between nodes 'n' and 'nn' during time period 't' assuming pipe diameter 'd' [MPa]")
        deltaPgas = Variable(m, "deltaPgas", type="positive", domain=[n, nn, t], description="Pressure drop 'gas-only' between nodes 'n' and 'nn' during time period 't' assuming pipe diameter 'd' [MPa]")
        deltaPliq = Variable(m, "deltaPliq", type="positive", domain=[n, nn, d, t], description="Pressure drop 'liquid-only' between nodes 'n' and 'nn' during time period 't' [MPa]")
        vel_liq = Variable(m, "vel_liq", type="positive", domain=[n, nn, d, t], description="Liquid velocity between nodes 'n' and 'nn' during time period 't' assuming pipe diameter 'd' [m/s]")

        accumulated_capacity = Variable(m, "accumulated_capacity", type="positive", domain=[pf, t, c], description="Total accumulated capacity at processing facility 'pf' of component 'c' during time period 't' [mscf per day]")
        pipe_cost = Variable(m, "pipe_cost", type="positive", domain=t, description="Total cost on pipeline installation during time period 't' [kUSD]")
        facility_cost = Variable(m, "facility_cost", type="positive", domain=t, description="Total cost on facility installation during time period 't' [kUSD]")
        total_cost = Variable(m, "total_cost", type="free", description="Total discounted cost [MMUSD]")

        # Call equations
        eqs = []
        eqs += mass_balances(m)
        eqs += costs_computation(m)
        eqs += capacity_constraints(m)
        eqs += pressure_bounds(m)
        eqs += liquid_only_pressure_drop(m)
        eqs += gas_only_pressure_drop(m)
        eqs += lm_correlation(m)

        model = Model(
            m,
            "Multiphase_network_design",
            equations=eqs,
            # equations=m.getEquations(),
            sense=Sense.MIN,
            problem=Problem.MIQCP,
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

        # Hotfix: Allow all connections
        arcs = []
        for i in self.m["i"].records["n"]:
            for j in self.m["j"].records["n"]:
                arcs.append((i, j))

        for j in self.m["j"].records["n"]:
            for pf in self.m["pf"].records["n"]:
                arcs.append((j, pf))


        self.m["arcs"].setRecords(arcs)
        self.m["tp"].setRecords([f"t{int(i)}" for i in list(self.m["st_time"].records["value"].unique())])
        self.m["dp_max"].setRecords(self.m["pwell"].records["value"].max() - self.m["pmin_pf"].records["value"].min())  # Set default max pressure drop
        self.m["maxPress"].setRecords(self.m["pwell"].records["value"].max())
        self.instance_maxFlow()
        self.compute_first_echelon_parameters()
        self.define_weymouth_constant()
        i = self.m["i"]
        self.m["Qinter"].fx[i, self.m["nn"], self.m["d"], self.m["t"], self.m["c"]].where[Ord(self.m["t"]) < self.m["st_time"][i]] = 0


    def instance_maxFlow(self):
        """
        Placeholder for now
        """
        maxFlow = self.m["maxFlow"]
        q_prod = self.m["Qprod"]
        i = self.m["i"]
        c = self.m["c"]
        t = self.m["t"]
        maxFlow[c,t] = Sum(i, q_prod[c, t])

    def compute_first_echelon_parameters(self):
        # TODO: Profile and improve algorithm
        mpa_to_psi = 145.038
        dp_records = []
        fix_press_records = []

        df_st_time = self.m["st_time"].records
        df_Qprod = self.m["Qprod"].records
        df_pwell = self.m["pwell"].records
        df_dist = self.m["dist"].records
        df_diam = self.m["diam"].records

        d_list = list(self.m["d"].records["uni"])
        t_list = list(self.m["tp"].records["t"])
        df_arcs = self.m["arcs"].records[["n","nn"]]
        arcs_list = [(row.n, row.nn) for row in df_arcs.itertuples(index=False)]
        first_echelon_arcs = [t for t in arcs_list if t[0].startswith("i")]
        for ((i,j), d, t) in product(first_echelon_arcs, d_list, t_list):
            dist = df_dist[(df_dist["n"] == i) & (df_dist["nn"] == j)].value.sum()
            st = int(df_st_time[df_st_time["i"] == i].value.iloc[0])
            t_adjusted = int(t[1:]) + 1 - st
            if dist > 0 and t_adjusted >= 1:
                Q_aux = df_Qprod[df_Qprod["t"] == f"t{t_adjusted}"]
                Qoil = Q_aux[Q_aux["c"] == "oil"].value.sum()
                Qgas = Q_aux[Q_aux["c"] == "gas"].value.sum()
                Qwater = Q_aux[Q_aux["c"] == "water"].value.sum()
                p_inlet = df_pwell[df_pwell["t"] == f"t{t_adjusted}"].value.sum()
                diam = df_diam[df_diam["d"] == d].value.sum()
                dp, _, _ = compute_multiphase_pressure_drop(Qoil=Qoil, Qgas=Qgas, Qwater=Qwater, p_inlet=p_inlet*mpa_to_psi, dist=dist, diam=diam)
                fix_press = p_inlet - dp
            else:
                dp = 0
                fix_press = self.m["maxPress"].records["value"].max()
            fix_press = max(fix_press, 0)  # Ensure non-negative
            fix_press_records.append((i,j,d,t,fix_press))
            dp_records.append((i,j,d,t,dp))

        self.m["dp_firstEch"].setRecords(dp_records)
        self.m["fixPress"].setRecords(fix_press_records)

    def define_weymouth_constant(self):
        df_dist = self.m["dist"].records
        df_diam = self.m["diam"].records
        df_arcs = self.m["arcs"].records[["n","nn"]]
        arcs_list = [(row.n, row.nn) for row in df_arcs.itertuples(index=False)]
        second_echelon_arcs = [t for t in arcs_list if t[0].startswith("j")]
        d_list = list(self.m["d"].records["uni"])

        kw_records = []
        for ((j,pf), d) in product(second_echelon_arcs, d_list):
            dist = df_dist[(df_dist["n"] == j) & (df_dist["nn"] == pf)].value.sum()
            diam = df_diam[df_diam["d"] == d].value.sum()
            wey_const = compute_weymouth_constant(dist=dist, diam=diam)
            kw_records.append((j, pf, d, wey_const))

        self.m["kw"].setRecords(kw_records)

    def solve(self):
        ...

    def predict(self, input_data):
        print(f"Predicting with {self.model_name} for input: {input_data}")
        return "prediction_result"