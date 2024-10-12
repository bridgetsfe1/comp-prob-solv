import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Read and extract csv data
data= pd.read_csv("U_Cv_values.csv")   
Ts = data.iloc[:, 1]  
Cv = data.iloc[:, 3]  

# Find the temperature at which Cv is maximum 
T_dissociation = Ts[np.argmax(Cv)]

# Plot Cv vs Temperature
plt.plot(Ts, Cv, label='$C_v(T)$')
plt.xlabel('Temperature ($K$)')
plt.ylabel('Heat Capacity Cv ($J/K$)')
plt.axvline(T_dissociation, color='orange', linestyle='dashdot', label=f'Disscociation at T={T_dissociation:.2f}K')
plt.title('Heat Capacity vs Temperature for Argon Dimer')
plt.legend()
plt.grid()
plt.savefig("Cv_vs_T")