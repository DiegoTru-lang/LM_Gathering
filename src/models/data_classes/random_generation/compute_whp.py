import numpy as np

def estimate_whp(oil, gas, water, alpha=0.02, beta=0.005, P_min=100, noise=True):
    """
    Estimate wellhead pressure from production values.
    
    Parameters:
        oil   : float, oil production (e.g., bbl/day)
        gas   : float, gas production (e.g., MSCF/day)
        water : float, water production (e.g., bbl/day)
        alpha : float, scaling factor for oil+water
        beta  : float, scaling factor for gas
        P_min : float, minimum wellhead pressure
        noise : bool, if True add small random noise
    
    Returns:
        whp : float, estimated wellhead pressure
    """

    pressure = P_min + alpha * (oil + water) + beta * gas
    if noise:
        pressure += np.random.normal(0, 5)  # ±5 psi random noise
    return max(pressure, 0)  # avoid negative pressures




if __name__ == "__main__":
    # Example: 500 bbl/d oil, 1000 Mcf/d gas, 200 bbl/d water
    p_wh = estimate_whp(oil=5767/2, gas=13304/2, water=11465/2)
    psi_to_mpa = 0.00689476
    print(f"Estimated wellhead pressure = {p_wh:.2f} psi")
    print(f"Estimated wellhead pressure = {p_wh*psi_to_mpa:.2f} MPa")