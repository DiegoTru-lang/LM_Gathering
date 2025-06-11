from enum import Enum

class ModelSetMapping(Enum):
    # Map model set/parameter names to DataClass keys
    I = ("i", "I")  # (model_set_name, data_key)
    J = ("j", "J")
    PF = ("pf", "PF")
    T = ("t", "T")
    C = ("c", "C")
    D = ("d", "D")
    S = ("s", "S")
    N = ("n", "N")
    ARCS = ("arcs", "ARCS")
    # Add more as needed
