import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k

T = 300
beta = 1 / (k * T)

def thermo(Z, T):
    """Calculates Internal Energy, Free Energy, and Entropy based on the partition function of a system
    Parameters:
        Z: Partition function
    Returns:
        U: float, Internal Energy
        F: float, free energy
        S: float, Entropy
    """
    # Internal energy (U)
    U = np.gradient(np.log(Z), 1 / (k * T)) * (-1)
    #Free Energy (F)
    F = -k * T * np.log(Z)
    #Entropy (S)
    S = np.gradient(F, T) * (-1)
    return U, F, S

def isolated_partition_function(energies, degeneracies, T):
    Z_iso = 0.0

    for E, g in zip(energies, degeneracies):
        Z_iso += g * np.exp(-E / (k * T))
    
    return Z_iso

Z_iso = isolated_partition_function(energies = [0]*14, degeneracies= [1]*14, T=300)
thermo_iso = thermo(Z_iso, T)


def soc_partition_function(energies, degeneracies, T):
    Z_soc = 0.0

    for E, g in zip(energies, degeneracies):
        Z_soc += g * np.exp(-E / (k * T))
    
    return Z_soc

Z_soc = soc_partition_function(energies = [0]*6 + [4.48609*(10**-20)]*8, degeneracies = [1]*14, T=300)
thermo_soc = thermo(Z_soc, T)

def cfs_partition_function(energies, degeneracies, T):
    Z_cfs = 0.0

    for E, g in zip(energies, degeneracies):
        Z_cfs += g * np.exp(-E / (k * T))
    
    return Z_cfs

Z_cfs = cfs_partition_function(energies = [0]*4 + [1.92261*(10**-20)]*2 + [4.00554*(10**-20)]*2 + [5.12697*(10**-20)]*4 + [7.37001*(10**-20)]*2, degeneracies = [1]*14, T=300)
thermo_cfs = thermo(Z_cfs, T)


temperatures = np.linspace(300, 2000, 100) 

U_isolated = []
F_isolated = []
S_isolated = []

U_soc = []
F_soc = []
S_soc = []

U_cfs = []
F_cfs =[]
S_cfs =[]

for T in temperatures:
    U_i, F_i, S_i = thermo(Z_iso, T)
    U_isolated.append(U_i)
    F_isolated.append(F_i)
    S_isolated.append(S_i)

    U_s, F_s, S_s = thermo(Z_soc, T)
    U_soc.append(U_s)
    F_soc.append(F_s)
    S_soc.append(S_s)

    U_c, F_c, S_c = thermo(Z_cfs, T)
    U_cfs.append(U_c)
    F_cfs.append(F_c)
    S_cfs.append(S_c)    


#Plotting Internal Energy
plt.figure(figsize=(8, 6))
plt.plot(temperatures, U_isolated, label='Isolated Ce', color='blue')
plt.plot(temperatures, U_soc, label='Ce with SOC', color='orange')
plt.plot(temperatures, U_cfs, label='Ce with SOC + CFS', color='green')
plt.title('Internal Energy vs Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Internal Energy (J)')
plt.legend()
plt.grid()
plt.savefig("Internal_Energy.png")  
plt.show()

# Plotting Free Energy
plt.figure(figsize=(8, 6))
plt.plot(temperatures, F_isolated, label='Isolated Ce', color='blue')
plt.plot(temperatures, F_soc, label='Ce with SOC', color='orange')
plt.plot(temperatures, F_cfs, label='Ce with SOC + CFS', color='green')
plt.title('Free Energy vs Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Free Energy (J)')
plt.legend()
plt.grid()
plt.savefig("Free_Energy.png")  
plt.show()

# Plotting Entropy
plt.figure(figsize=(8, 6))
plt.plot(temperatures, S_isolated, label='Isolated Ce', color='blue')
plt.plot(temperatures, S_soc, label='Ce with SOC', color='orange')
plt.plot(temperatures, S_cfs, label='Ce with SOC + CFS', color='green')
plt.title('Entropy vs Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Entropy (J/K)')
plt.legend()
plt.grid()
plt.savefig("Entropy.png") 
plt.show()


