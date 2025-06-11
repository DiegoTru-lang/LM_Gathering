from enum import Enum
from typing import Tuple
import pandas as pd

class ParameterLocation(Enum):
    #PARAM = (sheet_name, start_row, end_row, start_col, end_col)
    I = ("Sets", 1, 20, 0, 0) 
    J = ("Sets", 1, 20, 1, 1) 
    PF = ("Sets", 1, 20, 2, 2) 
    D = ("Sets", 1, 20, 3, 3) 
    S = ("Sets", 1, 20, 3, 3) 
    T = ("Sets", 1, 20, 3, 3) 
    IXLM_MAX = ("LM", 1, 1, 1, 1)  # single cell
    INIT_INTERVALS = ("LM", 2, 2, 1, 1)  
    PWELL = ("Pressure", 2, 12, 0, 1)  # range of cells
    HFFL = ("Fluid_dynamics", 1, 1, 1, 1)  
    RHO_LIQ = ("Fluid_dynamics", 2, 2, 1, 1)  
    KW = ("Fluid_dynamics", 3, 3, 1, 1)  

def read_parameter_from_excel(excel_path: str, location: Tuple[str, int, int, int, int]):
    sheet_name, start_row, end_row, start_col, end_col = location
    df = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)
    data = df.iloc[start_row:end_row+1, start_col:end_col+1]
    if data.size == 1:
        return data.values[0][0]
    elif data.shape[1] == 1:
        return list(data.values.flatten())
    return list(map(tuple, data.values))
