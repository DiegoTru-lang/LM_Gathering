from .data_classes import DataClass
from .data_classes.utils.model_mappings import ModelSetMapping
from gamspy import Container, Set, Sum, Parameter, Variable, Alias, Model, Sense, Ord, Options, Problem
from .equations import *
from .fluid_dynamics import compute_multiphase_pressure_drop, compute_weymouth_constant
from itertools import product
import sys
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from math import pi, sqrt
from math import pi, sqrt

class MultiphaseNetworkModel():
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
        
        selected_pipes = Set(m, "sel_pipes", domain=[n,nn,d], description="Pipeline connections that were selected in previous iterations of the algorithm")
        pw = Set(m, "pw", description="Segments for LM correalation piecewise linearization")

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
        maxFlow = Parameter(m, "maxFlow", domain=[c], description="Maximum flow of component 'c' [mscf per day]") # TODO: Compute as sum of production at all source nodes
        # pwell = Parameter(m, "pwell", domain=[i,t], description="Wellhead pressure per source node 'i' at time period 't' [MPa]")
        dp_firstEch = Parameter(m, "dp_firstEch", domain=[i,j,d,t], description="Pre-computed pressure drop at first echelon for connection (i,j), diameter 'd' at time period 't' [MPa]")
        dp_max = Parameter(m, "dp_max", description="Maximum allowable pressure drop between two nodes [MPa]")
        pwell = Parameter(m, "pwell", domain=t, description="Wellhead pressure at every node at time period 't' [MPa]")
        pmin_pf = Parameter(m, "pmin_pf", description="Minimum inlet pressure at processing facility [MPa]")
        # fixPress = Parameter(m, "fixPress", domain=[n, nn, t], description="Pre-computed max pressure at junction 'j' during time period 't' if connection (i, j) is installed [MPa]")
        fixPress = Parameter(m, "fixPress", domain=[i, j, d, t], description="Pre-computed max pressure at junction 'j' during time period 't' if connection (i, j) is installed with diameter 'd' [MPa]")
        maxPress = Parameter(m, "maxPress", description="Max pressure at node [MPa]")

        ixlm_ub = Parameter(m, "ixlm_ub", domain=[pw, j, pf, d, t], description="IXLM upper bound for LM interval 'pw' for connection (j,pf,d) at time 't'")
        ylp = Parameter(m, "ylp", domain=[pw, j, pf, d, t], description="YLM multiplier for liquid phase for LM interval 'pw' for connection (j,pf,d) at time 't'")
        allowed_int = Parameter(m, "allowed_int", domain=[j,pf,d,t], description="Order of maximum allowed interval 'pw' ") #i.e., if 3, then pw = 1,2,3 allowed; if 1, then only pw = 1 allowed

        # Variables

        # These variables must instantiate just for [n, nn] included in arcs set 
        y_pf = Variable(m, "y_pf", domain=[pf, s, t], type="binary", 
                        description="Equals 1 if a processing facility of size 's' is installed at node 'pf' during time period 't'")
        x_bar = Variable(m, "x_bar", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        x_pw = Variable(m, "x_pw", domain=[j, pf, d, t, pw], type="binary", description="Equals 1 if connection 'j' to 'pf' with diameter 'd' lies between ixlm interval of piecewise linearization 'pw' ")
        # xf = Variable(m, "xf", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' transports fluid from node 'n' to 'nn' at time period 't'")
        # xr = Variable(m, "xr", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' transports fluid from node 'nn' to 'n' at time period 't'")
        # xnf = Variable(m, "xnf", domain=[n, nn, d, t], type="binary", description="Equals 1 if no pipeline segment of diameter 'd' is installed between 'n' and 'nn' and P(n) <= P(nn) at time period 't'")
        # xnr = Variable(m, "xnr", domain=[n, nn, d, t], type="binary", description="Equals 1 if a pipeline segment of diameter 'd' between nodes 'n' and 'nn' is installed at time period 't'")
        
        q_inter = Variable(m, "Qinter", type="positive", domain=[n, nn, d, t, c], description="Flow of component 'c' through pipeline segment between nodes 'n' and 'nn' of diameter 'd' during time period 't' [mscf per day]")
        q_process = Variable(m, "Qprocess", type="positive", domain=[pf, t, c], description="Amount of component 'c' processed at facility 'pf' during time period 't' [mscf per day]")
        press = Variable(m, "press", type="positive", domain=[n, t], description="Pressure at node 'n' during time period 't' [MPa]")
        pressSQ = Variable(m, "pressSQ", type="positive", domain=[n, t], description="Squared pressure at node 'n' during time period 't' [MPa]")
        # pressGAS = Variable(m, "pressGAS", type="positive", domain=[n, t], description="Pressure at node 'n' during time period 't' assuming gas-only pressure drop [MPa]")
        pressGAS = Variable(m, "pressGAS", type="positive", domain=[pf, t], description="Pressure at node 'pf' during time period 't' assuming gas-only pressure drop [MPa]")
        deltaP = Variable(m, "deltaP", type="positive", domain=[n, nn, t], description="Pressure drop 'multiphase' between nodes 'n' and 'nn' during time period 't' [MPa]")
        deltaPgas = Variable(m, "deltaPgas", type="positive", domain=[n, nn, t], description="Pressure drop 'gas-only' between nodes 'n' and 'nn' during time period 't' assuming pipe diameter 'd' [MPa]")
        deltaPliq = Variable(m, "deltaPliq", type="positive", domain=[n, nn, d, t], description="Pressure drop 'liquid-only' between nodes 'n' and 'nn' during time period 't' [MPa]")
        deltaPliqYLM = Variable(m, "deltaPliqYLM", type="positive", domain=[n, nn, d, t, pw], description="Intermediate variable equal to pressure drop 'liquid-only' between nodes 'n' and 'nn' with diameter 'd' during time period 't' [MPa]")
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
        
        ## Initialize piecewise linearization
        j = self.m["j"]
        pf = self.m["pf"]
        t = self.m["t"]
        tp = self.m["tp"]
        d = self.m["d"]
        pw = self.m["pw"]
        ixlm_ub = self.m["ixlm_ub"]
        ylp = self.m["ylp"]
        allowed_int = self.m["allowed_int"]

        pw.setRecords([f"pw{p}" for p in range(1,11)])
        # ixlm_ub.setRecords([("pw1", 1)])
        # ylp.setRecords([("pw1", 1)])
        ixlm_ub["pw1", j, pf, d, t] = 3
        ylp["pw1", j, pf, d, t] = 1
        # ixlm_ub[pw, j, pf, d, tp].where[Ord(pw) > 1] = 5
        # ylp[pw, j, pf, d, tp].where[Ord(pw) > 1] = 5
        allowed_int[j, pf, d, t] = 1

    def instance_maxFlow(self):
        """
        Placeholder for now
        """
        maxFlow = self.m["maxFlow"]
        q_prod = self.m["Qprod"]
        i = self.m["i"]
        c = self.m["c"]
        t = self.m["t"]
        maxFlow[c] = Sum(i, q_prod[c, "t1"])

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
                dp, _, _, _, _ = compute_multiphase_pressure_drop(Qoil=Qoil, Qgas=Qgas, Qwater=Qwater, p_inlet=p_inlet*mpa_to_psi, dist=dist, diam=diam)
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

    def solve(self, solver: str = "gurobi", gap: float = 0.0001, max_time: float = 300.0):
        model = self.m.models["Multiphase_network_design"]
        summary = model.solve(solver=solver, options=Options(relative_optimality_gap=gap, time_limit = max_time))
        return summary
        # TODO: Validar si los cambios (valores en variables) se reflejan en self.m

    def plot_network(self):
        #TODO: Fix positions based on loc
        G = nx.DiGraph()

        source_nodes = self.m["i"].records["n"].tolist()
        junction_nodes = self.m["j"].records["n"].tolist()
        pf_nodes = self.m["pf"].records["n"].tolist()

        node_colors = []

        for _, row in self.m["n"].records.iterrows():
            G.add_node(row["n"], pos=(self.m["loc_x"][row["n"]], self.m["loc_y"][row["n"]]))
            if row["n"] in source_nodes:
                node_colors.append('lightblue')
            elif row["n"] in junction_nodes:
                node_colors.append('lightgreen')
            elif row["n"] in pf_nodes:
                node_colors.append('lightcoral')

        # Add edges based on installed pipelines
        sel_pipes = self.obtain_var_df("x_bar")
        sel_pipes = sel_pipes[sel_pipes["level"] > 0.5]
        
        for _, row in sel_pipes.iterrows():
            G.add_edge(row["n"], row["nn"], diameter=row["d"])

        pos = nx.get_node_attributes(G, 'pos')
        plt.figure(figsize=(10, 8))
        
        # Get edge widths based on diameter
        edge_widths = []
        for u, v, data in G.edges(data=True):
            try:
                width = float(data['diameter'])
            except Exception:
                width = 1.0
                edge_widths.append(width)

        nx.draw(
            G, pos, with_labels=True, node_size=500, node_color='lightblue',
            font_size=10, font_weight='bold', arrows=True, width=edge_widths
        )
        
        # Draw edge labels (diameters)
        edge_labels = {(u, v): f"d={data['diameter']}" for u, v, data in G.edges(data=True)}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

        plt.title("Pipeline Network")
        plt.show()

    def export_gdx(self, path: str = "toy_problem", run_name: str = "first_run"):
        model = self.m.models["Multiphase_network_design"]
        output_path = f"{path}/{run_name}"
        model.toGams(path=output_path)
        # TODO: export plot to same path

    def obtain_var_df(self, var_name: str, value_col: str = "level") -> pd.DataFrame:
        var_df = self.m[var_name].records
        idx = var_df.columns.get_loc(value_col)
        return var_df.iloc[:, :idx+1]
    
    def obtain_record_values(self, var_name: str, records: list[tuple], var_bool: bool = True, var_df: pd.DataFrame = None) -> dict[tuple, float | None]:
        value_column = "level" if var_bool else "value"
        if var_df is None:
            var_df = self.obtain_var_df(var_name, value_col=value_column)

        idx = var_df.columns.get_loc(value_column)
        key_cols = var_df.columns[:idx]

        results = {}
        for record in records:
            if len(record) != len(key_cols):
                raise ValueError(
                    f"Record {record} length ({len(record)}) does not match number of key columns ({len(key_cols)})."
                )

            mask = (var_df[key_cols] == record).all(axis=1)
            row = var_df.loc[mask, value_column]

            results[record] = row.iloc[0] if not row.empty else None

        return results
    
        # rec_df = pd.DataFrame(records, columns=key_cols)

        # # merge to find matches
        # merged = rec_df.merge(var_df, on=list(key_cols), how="left")

        # # build dictionary: tuple(record) -> level (or None if missing)
        # return {
        #     tuple(row[key_cols]): (row[value_column] if pd.notna(row[value_column]) else None)
        #     for _, row in merged.iterrows()
        # }

    def update_selected_pipes(self):
        selected_pipes = self.m["sel_pipes"]
        j = self.m["j"]

        sel_connect = self.obtain_var_df("x_bar")
        sel_connect = sel_connect[sel_connect["level"] > 0.5]

        for _, row in sel_connect.iterrows():
            if row["n"] not in j.records["n"].values:
                continue
            else:
                selected_pipes[row["n"], row["nn"], row["d"]] = True

        return sel_connect
    
    def compute_ixlm(self, qoil: float, qgas: float, qwater: float, p_inlet: float, dist: float, diam: float, kw_gas, n: float = 4.12):
        mpa_to_psi = 145.038
        m3_to_cf = 35.3147
        bbl_to_m3 = 0.158987
        in_to_m = 0.0254
        mile_to_km = 1.60934
        hffl = self.m["hffl"].records.value[0]
        rho_liq = self.m["rho_liq"].records.value[0]
        min_pressure = self.m["pmin_pf"].records.value[0]

        # dp, ixlm, ylp, dp_gas, dp_liq = compute_multiphase_pressure_drop(
        #     Qoil=qoil, Qgas=qgas, Qwater=qwater, p_inlet=p_inlet*mpa_to_psi, dist=dist, diam=diam
        # )
        
        ### Compute dpLiq
        vel_liq = ((qoil + qwater)*bbl_to_m3)/((pi*((diam*in_to_m)**2)/4)*24*3600)
        dp_liq = ((dist*mile_to_km)*(hffl/2)*(rho_liq*vel_liq*vel_liq)/(diam*in_to_m))/1e3

        ### Compute dpGas
        pressSQ = (p_inlet*1e3)**2
        aux = (pressSQ - ((qgas*1e3/m3_to_cf)*kw_gas)**2)
        if aux > 0:
            pressGAS = sqrt(aux)/1e3
        else:
            # pressGAS = min_pressure
            pressGAS = 0

        dp_gas = p_inlet - pressGAS

        ixlm = dp_gas/dp_liq
        ylp = ((ixlm**(1/n))+1)**n
        dp = dp_liq*ylp

        print(f"Inputs: qoil={qoil}, qgas={qgas}, qwater={qwater}, dist={dist}, diam={diam}, kw_gas={kw_gas}")
        print(f"Computed values: dp={dp:.2f}, ixlm={ixlm:.2f}, ylp={ylp:.2f}, dp_gas={dp_gas:.2f}, dp_liq={dp_liq:.2f}")
        print(f"p_inlet={p_inlet}, dp={dp:.2f}, press={p_inlet - dp}")

        max_incr = p_inlet - min_pressure
        ixlm_inf = dp_gas/(max_incr/ylp + dp_liq)
        ylp_inf = ((ixlm_inf**(1/n))+1)**n

        return (dp, ixlm, ixlm_inf, ylp, ylp_inf)

    def update_ixlm_intervals(self, last_connections: pd.DataFrame):
        i = self.m["i"]
        j = self.m["j"]
        ixlm_ub = self.m["ixlm_ub"]
        ixlm_df = self.obtain_var_df("ixlm_ub", value_col="value")
        ixlm_dict = {tuple(row[:-1]): row[-1] for row in ixlm_df.itertuples(index=False, name=None)}
        ylp = self.m["ylp"]
        ylp_df = self.obtain_var_df("ylp", value_col="value")
        ylp_dict = {tuple(row[:-1]): row[-1] for row in ylp_df.itertuples(index=False, name=None)}
        allowed_int = self.m["allowed_int"]
        allowed_int_df = self.obtain_var_df("allowed_int", value_col="value")
        allowed_int_dict = {tuple(row[:-1]): row[-1] for row in allowed_int_df.itertuples(index=False, name=None)}
        pw_list = self.m["pw"].records["uni"]
        c_list = self.m["c"].records["uni"]
        tp_list = self.m["tp"].records["t"]
        dist_df = self.obtain_var_df("dist", value_col="value")
        diams_df = self.obtain_var_df("diam", value_col="value")
        feasibility = True
        fixPress = self.m["fixPress"]
        fixPress_dict = {(row["i"], row["j"], row["d"], row["t"]): row["value"] for _, row in fixPress.records.iterrows()}
        init_pressure_dict = {(j, t): 10 for j in j.records["n"].values for t in tp_list}
        for _, row in last_connections.iterrows():
            if row["n"] in i.records["n"].values:
                for t in tp_list:
                    p_aux = fixPress_dict[(row["n"], row["nn"], row["d"], t)]
                    if p_aux < init_pressure_dict[(row["nn"], t)]:
                        init_pressure_dict[(row["nn"], t)] = p_aux
        
        for _, row in last_connections.iterrows():
            if not row["n"] in j.records["n"].values:
                continue
            j_aux, pf_aux, d_aux = (row["n"], row["nn"], row["d"])

            print(f"Processing node {j_aux} to {pf_aux} with diameter {d_aux}")

            t_aux = row["t"]
            dist = dist_df.loc[(dist_df['n'] == j_aux) & (dist_df['nn'] == pf_aux), "value"].values[0]
            diam = diams_df.loc[diams_df["d"] == d_aux, "value"].values[0]
            records = [(j_aux, pf_aux, d_aux, t, c) for t in tp_list for c in c_list]
            # TODO: Debería obtener los dataframes previo a las iteraciones
            press_values = self.obtain_record_values("press", [(j_aux, t) for t in tp_list])
            qinter_values = self.obtain_record_values("Qinter", records)
            
            for t in tp_list:
                max_int = allowed_int_dict[(j_aux, pf_aux, d_aux, t)]
                if int(t[1:]) >= int(t_aux[1:]):
                    qoil = qinter_values[(j_aux, pf_aux, d_aux, t, "oil")]
                    qgas = qinter_values[(j_aux, pf_aux, d_aux, t, "gas")]
                    qwater = qinter_values[(j_aux, pf_aux, d_aux, t, "water")]
                    print(f"  Time {t}: qoil={qoil}, qgas={qgas}, qwater={qwater}")
                    if (qoil is None) or (qoil == 0) or (qgas is None) or (qgas == 0) or (qwater is None) or (qwater == 0):
                        continue
                    # p_inlet = press_values[(j_aux, t)]
                    p_inlet = init_pressure_dict[(j_aux, t)]
                    kw_gas = self.obtain_record_values("kw", [(j_aux, pf_aux, d_aux)], var_bool=False)[(j_aux, pf_aux, d_aux)]
                    print("FLAG")
                    dp, ixlm_new, ixlm_inf, ylp_new, ylp_inf = self.compute_ixlm(qoil=qoil, qgas=qgas, qwater=qwater, p_inlet=p_inlet, dist=dist, diam=diam, kw_gas=kw_gas)
                    if p_inlet - dp < self.m["pmin_pf"].records.value[0]:
                        feasibility = False
                    if max_int > 1:
                        ixlm_values = self.obtain_record_values("ixlm_ub", [(pw, j_aux, pf_aux, d_aux, t) for pw in [f"pw{i+1}" for i in range(int(max_int) - 1)]], var_bool=False, var_df=ixlm_df)
                        ylp_values = self.obtain_record_values("ylp", [(pw, j_aux, pf_aux, d_aux, t) for pw in [f"pw{i+1}" for i in range(int(max_int) - 1)]], var_bool=False, var_df=ylp_df)
                        ixlm_list = [v for k, v in ixlm_values.items()]
                        ylp_list = [v for k, v in ylp_values.items()] #TODO: Verify it the lists are in the same order
                        # ixlm_list += [ixlm_inf, ixlm_new]
                        ixlm_list += [ixlm_new]
                        # ylp_list += [ylp_inf, ylp_new]
                        ylp_list += [ylp_new]
                    else:
                        ixlm_list = [ixlm_inf, ixlm_new]
                        ylp_list = [ylp_inf, ylp_new]
                    combined = sorted(zip(ixlm_list, ylp_list), reverse=False)
                    print(combined)
                    for i in range(len(combined)):
                        ixlm_dict[(f"pw{i+1}", j_aux, pf_aux, d_aux, t)] = float(combined[i][0])
                        # records_dict[ixlm_ub].append((f"pw{i+1}", j_aux, pf_aux, d_aux, t, float(combined[i][0]) if i < len(combined) else 3))
                        # ylp_dict[(f"pw{i+1}", j_aux, pf_aux, d_aux, t)] = float(combined[i-1][1]) if i > 0 else 1
                        ylp_dict[(f"pw{i+1}", j_aux, pf_aux, d_aux, t)] = float(combined[i][1])
                        # records_dict[ylp].append((f"pw{i+1}", j_aux, pf_aux, d_aux, t, float(combined[i-1][1]) if i > 0 else 1))
                    allowed_int_dict[(j_aux, pf_aux, d_aux, t)] = len(combined) + 1
                    # records_dict[allowed_int].append((j_aux, pf_aux, d_aux, t, len(combined) + 1))
        records_dict = {
            ixlm_ub: [(k[0], k[1], k[2], k[3], k[4], v) for k, v in ixlm_dict.items()],
            ylp: [(k[0], k[1], k[2], k[3], k[4], v) for k, v in ylp_dict.items()],
            allowed_int: [(k[0], k[1], k[2], k[3], v) for k, v in allowed_int_dict.items()]
        }
        return (feasibility, records_dict)

    def solution_algorithm(self, scenario_name: str, max_iterations: int = 5, solver_opts: dict = {"solver": "gurobi", "gap": 0.0001, "max_time": 300.0}):
        it = 0
        feasibility = False
        while not feasibility and it < max_iterations:
            print(f"--- Iteration {it+1} ---")
            self.solve(solver=solver_opts.get("solver", "gurobi"), gap=solver_opts.get("gap", 0.0001), max_time=solver_opts.get("max_time", 300.0))
            # self.export_gdx(path=scenario_name, run_name=f"iteration_{it+1}")
            sel_connect = self.update_selected_pipes()
            feasibility, records_dict = self.update_ixlm_intervals(sel_connect)
            self.m.setRecords(records_dict)
            self.export_gdx(path=scenario_name, run_name=f"iteration_{it+1}")
            it += 1

        return (feasibility, it)

    def predict(self, input_data):
        print(f"Predicting with {self.model_name} for input: {input_data}")
        return "prediction_result"