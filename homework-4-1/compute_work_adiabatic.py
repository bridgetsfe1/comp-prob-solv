from scipy.integrate import trapezoid
import numpy as np
import pandas as pd

n = 1
R = 8.314
T = 300
V_init = 0.1
gamma = 1.4

def work_adi(V_init, V_final, T, gamma):
    """Calculates work in an adiabatic expansion.
    Parameters:
        V_init: intial volume
        V_final: final volume
        T: temperature
        gamma: heat capacity ratio
    Returns:
        float: work done adiabatically
    """
    P_init = (n*R*T)/ V_init
    C = P_init * (V_init ** gamma)
    V = np.linspace(V_init, V_final, 100)
    P = C / (V**gamma)
    work = -trapezoid(P, V)
    return work

Vf_values = np.linspace(V_init, 3 * V_init, num=50)
adiabatic_work = [work_adi(V_init, Vf, T, gamma) for Vf in Vf_values]

rows = {
    'Final Volume (Vf) [m³]': Vf_values,
    'Adiabatic Work (J)': adiabatic_work,
}
df = pd.DataFrame(rows)

df.to_csv('adi_work_vs_final_volume.csv')


