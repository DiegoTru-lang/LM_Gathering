from enum import Enum
from typing import Tuple
import pandas as pd

# Enum to specify parameter locations in the Excel file
class ParameterLocation(Enum):
    # Example: PARAM1 = (sheet_name, start_row, end_row, start_col, end_col)
    PARAM1 = ("Sheet1", 1, 1, 0, 0)  # single cell
    PARAM2 = ("Sheet1", 2, 5, 1, 1)  # range of cells
    # Add more parameters as needed

def read_parameter_from_excel(excel_path: str, location: Tuple[str, int, int, int, int]):
    sheet_name, start_row, end_row, start_col, end_col = location
    df = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)
    data = df.iloc[start_row:end_row+1, start_col:end_col+1]
    if data.size == 1:
        return data.values[0][0]
    return data.values.tolist()
