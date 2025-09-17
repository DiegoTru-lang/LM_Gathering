from .NetworkModel import MultiphaseNetworkModel
from .data_classes import DataClass
from .data_classes.utils.model_mappings import ModelSetMapping
from gamspy import Container, Set, Sum, Parameter, Variable, Alias, Model, Sense, Ord, Options, Problem
from .equations import *
from .fluid_dynamics import compute_multiphase_pressure_drop, compute_weymouth_constant
import sys

class GatheringModel(MultiphaseNetworkModel):
    def __init__(self, model_name: str, data: DataClass = None, fixed_diameters: bool = True):
        super().__init__(model_name, data)
        self.fixed_diameters = fixed_diameters
        if fixed_diameters:
            pass
        # Two options for the model:
        # 1. With fixed diameters
        # 2. With diameter selection:
        #   a. For j -> pf only
        #   b. For every i -> j and j -> pf

    def _construct_model(self) -> Container:
        m = Container(output=sys.stdout)

        # Sets

        n = Set(m, "n", description="Nodes in the network")
        nn = Alias(m, name="nn", alias_with=n)
        i = Set(m, "i", domain=n, description="Source nodes / multi-phase")
        j = Set(m, "j", domain=n, description="Junction nodes / multi-phase")
        pf = Set(m, "pf", domain=n, description="Processing facility")

        arcs = Set(m, "arcs", domain=[n,nn], description="Allowed connections")
        
        t = Set(m, "t", description="Time periods")
        tt = Alias(m, name="tt", alias_with=t)
        tp = Set(m, "tp", domain=t,description="Time periods where there is a production peak")
        c = Set(m, "c", records=["oil", "gas", "water"], description="Flow components")
        d = Set(m, "d", description="Pipeline diameter options")
        s = Set(m, "s", description="Facility sizes")
        
        pw = Set(m, "pw", description="Segments for LM correalation piecewise linearization")

        # Parameters

        # q_prod = Parameter(m, "Qprod", domain=[c,t], description="Oil production at each node at 't' periods after start-time [BBL per day]")
        q_prod = Parameter(m, "Qprod", domain=[i,t,c], description="Production at source node 'i' of component 'c' during time period 't' [mscf per day]")
        st_time = Parameter(m, "st_time", domain=i, description="Production start time of source node 'i' ") #TODO: Compute based on q_prod
        diam = Parameter(m, "diam", domain=d, description="Diameter of pipeline option 'd' [inches]")
        diameter_cost = Parameter(m, "diameter_cost", domain=d, description="Pipeline cost per km for diameter size 'd' [kUSD/mile]")
        
        ir = Parameter(m, "ir", description="Interest rate [%]")

        loc_x = Parameter(m, "loc_x", domain=n, description="Position of node 'n' in x-axis [km]")
        loc_y = Parameter(m, "loc_y", domain=n, description="Position of node 'n' in y-axis [km]")
        dist = Parameter(m, "dist", domain=[n,nn], description="Distance between nodes'n' to 'nn' [mile]")

        hffl = Parameter(m, "hffl", records=0.025, description="Hydraulic friction factor for liquid phase")
        rho_liq = Parameter(m, "rho_liq", records=850, description="Density of liquid phase mixture")
        kw = Parameter(m, "kw", domain=[n,nn,d], description="Weymouth constants and parameters synthesized")
        maxFlow = Parameter(m, "maxFlow", domain=[c,t], description="Maximum flow of component 'c' during time period 't' [mscf per day]") # TODO: Compute as sum of production at all source nodes
        pwell = Parameter(m, "pwell", domain=[i,t], description="Wellhead pressure per source node 'i' at time period 't' [MPa]")
        dp_firstEch = Parameter(m, "dp_firstEch", domain=[i,j,d,t], description="Pre-computed pressure drop at first echelon for connection (i,j), diameter 'd' at time period 't' [MPa]")
        dp_max = Parameter(m, "dp_max", description="Maximum allowable pressure drop between two nodes [MPa]")
        pmin_pf = Parameter(m, "pmin_pf", description="Minimum inlet pressure at processing facility [MPa]")
        fixPress = Parameter(m, "fixPress", domain=[i, j, d, t], description="Pre-computed max pressure at junction 'j' during time period 't' if connection (i, j) is installed with diameter 'd' [MPa]")
        maxPress = Parameter(m, "maxPress", description="Max pressure at node [MPa]")

        ixlm_ub = Parameter(m, "ixlm_ub", domain=[pw, j, pf, d, t], description="IXLM upper bound for LM interval 'pw' for connection (j,pf,d) at time 't'")
        ylp = Parameter(m, "ylp", domain=[pw, j, pf, d, t], description="YLM multiplier for liquid phase for LM interval 'pw' for connection (j,pf,d) at time 't'")
        allowed_int = Parameter(m, "allowed_int", domain=[j,pf,d,t], description="Order of maximum allowed interval 'pw' ") #i.e., if 3, then pw = 1,2,3 allowed; if 1, then only pw = 1 allowed

        # Variables

        # These variables must instantiate just for [n, nn] included in arcs set 
        x_bar = Variable(m, "x_bar", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        x_pw = Variable(m, "x_pw", domain=[j, pf, d, t, pw], type="binary", description="Equals 1 if connection 'j' to 'pf' with diameter 'd' lies between ixlm interval of piecewise linearization 'pw' ")

        q_inter = Variable(m, "Qinter", type="positive", domain=[n, nn, d, t, c], description="Flow of component 'c' through pipeline segment between nodes 'n' and 'nn' of diameter 'd' during time period 't' [mscf per day]")
        q_interSQ = Variable(m, "QinterSQ", type="positive", domain=[n, nn, d, t, c], description="Squared flow of 'gas' through pipeline segment between nodes 'n' and 'nn' of diameter 'd' during time period 't' [(mscf per day)**2]")
        press = Variable(m, "press", type="positive", domain=[n, t], description="Pressure at node 'n' during time period 't' [MPa]")
        pressSQ = Variable(m, "pressSQ", type="positive", domain=[n, t], description="Squared pressure at node 'n' during time period 't' [MPa]")
        pressGAS = Variable(m, "pressGAS", type="positive", domain=[pf, t], description="Pressure at node 'pf' during time period 't' assuming gas-only pressure drop [MPa]")
        deltaP = Variable(m, "deltaP", type="positive", domain=[n, nn, t], description="Pressure drop 'multiphase' between nodes 'n' and 'nn' during time period 't' [MPa]")
        deltaPgas = Variable(m, "deltaPgas", type="positive", domain=[n, nn, t], description="Pressure drop 'gas-only' between nodes 'n' and 'nn' during time period 't' assuming pipe diameter 'd' [MPa]")
        deltaPliq = Variable(m, "deltaPliq", type="positive", domain=[n, nn, d, t], description="Pressure drop 'liquid-only' between nodes 'n' and 'nn' during time period 't' [MPa]")
        deltaPliqYLM = Variable(m, "deltaPliqYLM", type="positive", domain=[n, nn, d, t, pw], description="Intermediate variable equal to pressure drop 'liquid-only' between nodes 'n' and 'nn' with diameter 'd' during time period 't' [MPa]")
        vel_liq = Variable(m, "vel_liq", type="positive", domain=[n, nn, d, t], description="Liquid velocity between nodes 'n' and 'nn' during time period 't' assuming pipe diameter 'd' [m/s]")

        q_processed = Variable(m, "Qprocessed", type="positive", domain=[t,c], description="Total production of component 'c' during time period 't' [bbl/mscf per day]")
        total_processed = Variable(m, "total_processed", type="positive", description="Total production of component 'c' during time period 't' [bbl/mscf per day]")
        
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
            objective=total_processed,
            # options=Options(mip="gurobi"),
        )

        return m