import numpy as np
import pandas as pd


import numpy as np

def generate_random_wells_dirichlet(
    base_proportions, base_oil_profile, n_wells,
    kappa=30, random_state=None,
    noise=False, noise_level=0.05
):
    """
    Generate synthetic well production profiles with STATIC random oil/gas/water ratios.
    
    Parameters
    ----------
    base_proportions : list or array of length 3
        Base proportions [oil, gas, water] that define the mean composition.
    base_oil_profile : list or array
        Production profile of oil for the time horizon (e.g. 10 years).
    n_wells : int
        Number of random wells to generate.
    kappa : float, optional
        Concentration parameter for Dirichlet distribution (higher = less variance).
    random_state : int, optional
        Seed for reproducibility.
    noise : bool, optional
        If True, adds random multiplicative noise to production values.
    noise_level : float, optional
        Standard deviation of Gaussian noise (in relative terms, e.g. 0.05 = 5%).
    
    Returns
    -------
    dict
        Dictionary with keys (t, well, component) → production value.
    """
    rng = np.random.default_rng(random_state)
    p_mean = np.array(base_proportions)
    
    # Total production = oil / base oil fraction
    base_oil_profile = np.array(base_oil_profile)
    total_profile = base_oil_profile / p_mean[0]
    
    # Dirichlet sampling
    alpha = p_mean * kappa
    compositions = rng.dirichlet(alpha, size=n_wells)
    
    # Build dictionary
    result = {}
    components = ["oil", "gas", "water"]
    time_horizon = len(base_oil_profile)
    
    for w in range(n_wells):
        for t in range(time_horizon):
            for i, comp in enumerate(components):
                value = compositions[w, i] * total_profile[t]
                if noise:
                    # Multiplicative noise ~ N(1, noise_level)
                    value *= rng.normal(loc=1.0, scale=noise_level)
                    value = max(value, 0)  # prevent negatives
                result[(t+1, f"i{w+1}", comp)] = value
    
    return result


if __name__ == "__main__":
    base_p = [0.15, 0.3, 0.55]
    oil_profile = [6500, 6100, 5200, 4400, 4000, 3700, 3400, 3200, 3000, 2800]

    res = generate_random_wells_dirichlet(
        base_p, oil_profile, n_wells=3,
        random_state=42, noise=True, noise_level=0.1
    )

    df = pd.Series(res).unstack(level=2)
    print(df.head(12))
