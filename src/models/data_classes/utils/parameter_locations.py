from enum import Enum
from typing import Tuple
import pandas as pd

class ParameterLocation(Enum):
    #PARAM = (sheet_name, start_row, end_row, start_col, end_col)
    I = (("Sets", 1, 20, 0, 0), False)
    J = (("Sets", 1, 20, 1, 1), False)
    PF = (("Sets", 1, 20, 2, 2), False)
    D = (("Sets", 1, 20, 3, 3), False)
    S = (("Sets", 1, 20, 4, 4), False)
    T = (("Sets", 1, 200, 5, 5), False)
    C = (("Sets", 1, 20, 6, 6), False)
    A = (("Sets", 1, 200, 8, 9), False)
    LOC_X = (("Location", 1, 50, 0, 1), False)
    LOC_Y = (("Location", 1, 50, 3, 4), False)
    DIST_IJ = (("Distances", 0, 12, 0, 12), "dist_ij")
    DIST_JK = (("Distances", 0, 12, 15, 18), "dist_jk")
    PROD = (("Production", 0, 200, 0, 3), "Qprod")
    START_TIME = (("Production", 1, 50, 8, 9), False)
    FAC_CAP = (("Capacity", 1, 20, 0, 2), False)
    DIAM_SIZES = (("Capacity", 1, 20, 7, 8), False)
    FAC_COST = (("Capital expenditure", 1, 20, 0, 1), False)
    DIAM_COST = (("Capital expenditure", 1, 20, 3, 4), False)
    IR = (("Capital expenditure", 1, 1, 7, 7), False)
    IXLM_MAX = (("LM", 1, 1, 1, 1), False)  # single cell
    INIT_INTERVALS = (("LM", 2, 2, 1, 1), False)
    PWELL = (("Pressure", 2, 200, 0, 1), False)  # range of cells
    P_IN = (("Pressure", 1, 1, 4, 4), False)
    P_OUT = (("Pressure", 2, 2, 4, 4), False)
    HFFL = (("Fluid_dynamics", 1, 1, 1, 1), False)
    RHO_LIQ = (("Fluid_dynamics", 2, 2, 1, 1), False) 
    KW = (("Fluid_dynamics", 3, 3, 1, 1), False)  
    GOR = (("Fluid_dynamics", 4, 4, 1, 1), False)  
    WOR = (("Fluid_dynamics", 5, 5, 1, 1), False)  

def read_parameter_from_excel(excel_path: str, location: Tuple[str, int, int, int, int], pivot_table: str = None):
    sheet_name, start_row, end_row, start_col, end_col = location
    if pivot_table:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
        data = df.iloc[start_row:end_row+1, start_col:end_col+1]
        if pivot_table == "Qprod":
            df_melted = data.melt(id_vars=["Time period"], 
                                    value_vars=["Oil production [BBL/day]", "Gas production [mscf/day]", "Water production [BBL/day]"],
                                    var_name="Component",
                                    value_name="Value")
            df_melted["Component"] = df_melted["Component"].map({
                                    "Oil production [BBL/day]": "oil",
                                    "Gas production [mscf/day]": "gas",
                                    "Water production [BBL/day]": "water"})
            df_melted = df_melted[["Component", "Time period", "Value"]]
        else:
            df_melted = data.melt(id_vars=[pivot_table], var_name="j", value_name="value")
            df_melted["value"] = df_melted["value"].fillna(0)
            df_melted = df_melted.dropna(subset=[pivot_table, "j"])
        return list(df_melted.itertuples(index=False, name=None))
    df = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)
    data = df.iloc[start_row:end_row+1, start_col:end_col+1]
    data = data.dropna()
    if data.size == 1:
        return data.values[0][0]
    elif data.shape[1] == 1:
        return list(data.values.flatten())
    return list(map(tuple, data.values))
