import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k
import pandas as pd

temperatures = np.linspace(300, 2000, 100) 
beta = 1 / (k * temperatures)

def thermo(Z, T):
    """Calculates Internal Energy, Free Energy, and Entropy based on the partition function of a system
    Parameters:
        Z: (nparray) Partition function
        T: (nparray) temperatures
    Returns:
        nparray: Internal Energy
        nparray: Free energy
        nparray: Entropy
    """
    # Internal energy (U)
    U = -np.gradient(np.log(Z), 1 / (k * T))
    # Free Energy (F)
    F = -k * T * np.log(Z)
    # Entropy (S)
    S = -np.gradient(F, T)
    return U, F, S

def isolated_partition_function(energies, degeneracies, T):
    """Calculates isolated partition function.
    Parameters: 
        energies (list): List of energy values
        degeneracies (list): List of degeneracies (equal value of energies)
        T (nparray): array of temepratures over which the partition fn is calculated. 
    Returns:
        nparray: Isolated partition function.
    """
    Z_iso = 0.0

    for E, g in zip(energies, degeneracies):
        Z_iso += g * np.exp(-E / (k * T))
    
    return Z_iso

Z_iso = isolated_partition_function(energies = [0]*14, degeneracies= [1]*14, T=temperatures)
U_i, F_i, S_i = thermo(Z_iso, temperatures)


def soc_partition_function(energies, degeneracies, T):
    """Calculates SOC partition function.
    Parameters: 
        energies (list): List of energy values
        degeneracies (list): List of degeneracies (equal value of energies)
        T (nparray): array of temepratures over which the partition fn is calculated. 
    Returns:
        nparray: SOC partition function.
    """
    Z_soc = 0.0

    for E, g in zip(energies, degeneracies):
        Z_soc += g * np.exp(-E / (k * T))
    
    return Z_soc

Z_soc = soc_partition_function(energies = [0]*6 + [4.48609*(10**-20)]*8, degeneracies = [1]*14, T=temperatures)
U_soc, F_soc, S_soc = thermo(Z_soc, temperatures)

def cfs_partition_function(energies, degeneracies, T):
    """Calculates SOC+CFS partition function.
    Parameters: 
        energies (list): List of energy values
        degeneracies (list): List of degeneracies (equal value of energies)
        T (nparray): array of temepratures over which the partition fn is calculated. 
    Returns:
        nparray: SOC+CFS partition function.
    """
    Z_cfs = 0.0

    for E, g in zip(energies, degeneracies):
        Z_cfs += g * np.exp(-E / (k * T))
    
    return Z_cfs

Z_cfs = cfs_partition_function(energies = [0]*4 + [1.92261*(10**-20)]*2 + [4.00554*(10**-20)]*2 + [5.12697*(10**-20)]*4 + [7.37001*(10**-20)]*2, degeneracies = [1]*14, T=temperatures)
U_cfs, F_cfs, S_cfs = thermo(Z_cfs, temperatures)

#All these functions do the same, could make one function for calculating the partition function and have different energy inputs (followed instructions to make 3 separate fns, but a way to improve efficiency!)

temperatures = np.linspace(300, 2000, 100)  

all_cases = {
    "U (iso)": U_i,
    "U (SOC)": U_soc,
    "U (SOC+CFS)":U_cfs,
    "F (iso)": F_i, 
    "F (SOC)": F_soc, 
    "F (SOC+CFS)": F_cfs,
    "S (iso)": S_i, 
    "S (SOC)": S_soc, 
    "S (SOC+CFS)":S_cfs
}

df = pd.DataFrame(all_cases, index=temperatures)
df.to_csv("thermo_prop_ce.csv")

#Plotting Internal Energy
plt.plot(temperatures, U_i, label='Isolated Ce', color='blue')
plt.plot(temperatures, U_soc, label='Ce with SOC', color='pink')
plt.plot(temperatures, U_cfs, label='Ce with SOC + CFS', color='purple')
plt.title('Internal Energy vs Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Internal Energy (J)')
plt.legend()
plt.grid()
plt.savefig("Internal_Energy.png")  
plt.show()
plt.close()

# Plotting Free Energy
plt.plot(temperatures, F_i, label='Isolated Ce', color='blue')
plt.plot(temperatures, F_soc, label='Ce with SOC', color='pink')
plt.plot(temperatures, F_cfs, label='Ce with SOC + CFS', color='purple')
plt.title('Free Energy vs Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Free Energy (J)')
plt.legend()
plt.grid()
plt.savefig("Free_Energy.png")  
plt.show()
plt.close()

# Plotting Entropy
plt.plot(temperatures, S_i, label='Isolated Ce', color='blue')
plt.plot(temperatures, S_soc, label='Ce with SOC', color='pink')
plt.plot(temperatures, S_cfs, label='Ce with SOC + CFS', color='purple')
plt.title('Entropy vs Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Entropy (J/K)')
plt.legend()
plt.grid()
plt.savefig("Entropy.png") 
plt.show()
plt.close()


