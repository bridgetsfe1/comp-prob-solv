import numpy as np
from scipy.constants import k, h
from scipy.integrate import trapezoid
import pandas as pd
import math

# Given Parameters (in Joules and meters)
epsilon = 1.650242e-21  
sigma = 3.4e-10  
V = 1000 * 1e-30  
T_min = 10
T_max = 1000
temps = np.linspace(T_min, T_max, 100) 
m = 39.948 * 1.6605e-27  # Mass of Argon in kg 
L_max = np.cbrt(V)

#Issue with thermal wavelength being super large
def thermal_wavelength(T):
    """
    Parameters:
        T (nparray): temperature 
    Returns:
        nparray: thermal wavelength
    """
    return math.sqrt((h ** 2) / (2 * np.pi * m * k * T)) 

#Copied from in-class code
def lennard_jones(r, epsilon, sigma):
    """
    Calculate the Lennard-Jones potential energy between two particles.
    Parameters:
        r (float): Distance between two particles.
        epsilon (float): Depth of the potential well (eV).
        sigma (float): Finite distance at which the inter-particle potential is zero (Angstrom).
    Returns:
        float: Potential energy (eV).
    """
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

# Partition function for two LJ particles in a cubic box
def partition_function(T):
    lamb = thermal_wavelength(T) 

    coeff = (4 * np.pi) / ((h ** 6) * (lamb** 6))  

    r_min = 0.001 * sigma  
    r_max = L_max
    
    # Discretize the distance space
    r_values = np.linspace(r_min, r_max, 1000)
    
    # Compute the Boltzmann factor for each relative distance r
    integrand_values = np.array([np.exp(-lennard_jones(r, epsilon, sigma) / (k * T)) * r**2 for r in r_values])
    

    # Perform the integration using the trapezoidal rule
    Z = trapezoid(integrand_values, r_values)
    
    return coeff * Z  

# Compute partition function for all temperatures
partition_values = []
for T in temps:
    partition_values.append(partition_function(T))

df = pd.DataFrame({"Temperature (K)": temps, "Partition Function": partition_values})
df.to_csv("partition_function_values.csv")




