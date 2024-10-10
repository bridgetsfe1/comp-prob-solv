from scipy.integrate import trapezoid
import numpy as np
import pandas as pd

n = 1
R = 8.314
T = 300
V_init = 0.1
gamma = 1.4

def work_iso(V_init, V_final, T):
    """Calculates work in an isothermal expansion
    Parameters:
        V_init: intial volume
        V_final: final volume
        T: constant temperature
    Returns:
        float: work done isothermally
    """
    V = np.linspace(V_init, V_final, 100)
    integrand = (n * R * T) / V
    work = -trapezoid(integrand, V)
    return work

Vf_values = np.linspace(V_init, 3 * V_init, num=50)
isothermal_work = [work_iso(V_init, Vf, T) for Vf in Vf_values]

rows = {
    'Final Volume (Vf) [m³]': Vf_values,
    'Isothermal Work (J)': isothermal_work,
}
df = pd.DataFrame(rows)

df.to_csv('iso_work_vs_final_volume.csv')

