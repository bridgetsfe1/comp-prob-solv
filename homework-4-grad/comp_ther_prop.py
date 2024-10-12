import numpy as np
from scipy.constants import k
import pandas as pd
import sys 
from comp_part_func import partition_function


#Given Parameters
epsilon = 1.650242e-21  
sigma = 3.4e-10  
V = 1000 * 1e-30  
temps = np.linspace(start = 10, stop = 1000, num= 10000) 
m = 39.948 * 1.6605e-27  # Mass of Argon in kg 
L_max = np.cbrt(V)


# Calculate internal energy U and heat capacity Cv
def U_and_Cv(temps):
    """"
    Computes internal energy and heat capacity using parition functions.
    Parameters:
        temps(array): Array of Temperatures
    Returns:
        float: U: internal energy
        float: Cv: heat capacity
    """
    #Create list to append partitions
    Z_values = []    
    
    # Calculate partition function and its logarithm for all temperatures
    for T in temps:
        Z_T = partition_function(T)
        Z_values.append(Z_T)
    
    Z = np.array(Z_values)
    
    beta = 1 / (k * temps)
    
    # Calculate internal energy U = -d(ln Z)/d beta
    U_values = -np.gradient(np.log(Z), beta)
    
    # Calculate heat capacity Cv = dU/dT
    U = np.array(U_values)
    Cv = np.gradient(U, temps)
    
    return U, Cv

# Calculate partition function, internal energy, and heat capacity
U, Cv = U_and_Cv(temps)

#Create data for csv file
df = pd.DataFrame({
    'Temperature (K)': temps,
    'Internal Energy (U)': U,
    'Heat Capacity (Cv)': Cv})
df.to_csv('U_Cv_values.csv')