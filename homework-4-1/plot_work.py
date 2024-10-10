from compute_work_isothermal import work_iso
from compute_work_adiabatic import work_adi
import matplotlib.pyplot as plt
import numpy as np

#Defining known values
n = 1
R = 8.314
T = 300
V_init = 0.1
gamma = 1.4
Vf_values = np.linspace(V_init, 3 * V_init, num=50)

isothermal_work = [work_iso(V_init, Vf, T) for Vf in Vf_values]
adiabatic_work = [work_adi(V_init, Vf, T, gamma) for Vf in Vf_values]

plt.plot(Vf_values, isothermal_work, label='Isothermal Work', color='blue')
plt.plot(Vf_values, adiabatic_work, label='Adiabatic Work', color='red')
plt.title('Work Done in Isothermal and Adiabatic Processes')
plt.xlabel('Final Volume ($m^3$)')
plt.ylabel('Work Done ($J$)')
plt.legend()
plt.grid()
plt.xlim(V_init, 3 * V_init)
plt.ylim(-3000, 0)
plt.savefig("Work_Figure")