import numpy as np
import pandas as pd

def generate_random_wells_logistic_normal(
    base_proportions, base_oil_profile, n_wells,
    phi=0.8, sigma=0.1, random_state=None
):
    """
    Generate synthetic well production profiles with time-varying proportions
    using a logistic-normal AR(1) process.
    
    Parameters
    ----------
    base_proportions : list of 3 floats
        Initial mean proportions [oil, gas, water].
    base_oil_profile : list or array
        Oil production profile for the time horizon.
    n_wells : int
        Number of wells to simulate.
    phi : float, optional
        AR(1) persistence parameter (0 < phi < 1).
    sigma : float, optional
        Std. dev. of the innovations in logit space.
    random_state : int, optional
        Seed for reproducibility.
    
    Returns
    -------
    dict
        Dictionary with keys (t, well, component) -> production value.
    """
    rng = np.random.default_rng(random_state)
    base_oil_profile = np.array(base_oil_profile)
    T = len(base_oil_profile)
    
    # Compute total profile from oil / oil fraction
    p0 = np.array(base_proportions)
    total_profile = base_oil_profile / p0[0]
    
    # Reference mean in logit space
    logit = lambda p: np.log(p[:-1] / p[-1])
    inv_logit = lambda z: np.append(
        np.exp(z) / (1 + np.sum(np.exp(z))),
        1.0 / (1 + np.sum(np.exp(z)))
    )
    mu = logit(p0)
    
    results = {}
    comps = ["oil", "gas", "water"]
    
    for w in range(n_wells):
        z_t = mu.copy()
        for t in range(T):
            # AR(1) step
            eps = rng.normal(0, sigma, size=len(z_t))
            z_t = mu + phi * (z_t - mu) + eps
            
            # Back-transform
            p_t = inv_logit(z_t)
            
            # Save production values
            for i, comp in enumerate(comps):
                results[(t+1, f"Well{w+1}", comp)] = p_t[i] * total_profile[t]
    
    return results

if __name__ == "__main__":
    # Example usage
    base_p = [0.15, 0.3, 0.55]
    oil_profile = [6500, 6100, 5200, 4400, 4000, 3700, 3400, 3200, 3000, 2800]

    res = generate_random_wells_logistic_normal(
        base_p, oil_profile, n_wells=3, phi=0.15, sigma=.15, random_state=42
    )

    df = pd.Series(res).unstack(level=2)  # easier to read
    print(df.head(15))
