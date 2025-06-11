from .parameter_locations import ParameterLocation, read_parameter_from_excel

def read_input(excel_filename: str = "parameters.xlsx"):
    """
    Reads all parameters defined in ParameterLocation from the Excel file.
    Returns a dictionary of parameter values.
    """
    import os
    params = {}
    excel_path = os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", excel_filename)
    for param in ParameterLocation:
        params[param.name] = read_parameter_from_excel(excel_path, param.value)
    return params