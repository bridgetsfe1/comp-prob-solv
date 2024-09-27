import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

#For sigma: 3.4 Angstroms = 3.4E-10 meters
sigma = 3.4
k_b = 8.617E-5      #eV/k
N_a = 6.022E23      #Avogadros number
#1. Hard Sphere Potential
def hard_sphere(r, sigma=3.4):
    """Calculates the hard sphere potential. 
    Parameters:
        r: distance between the atoms. 
        sigma: diameter of the hard sphere. 
    Returns:
        float: Potential energy V(r)
    """
    if r < sigma:
        return 1000  # Approximate infinity as 1000
    else:
        return 0  # Zero potential for distances greater than sigma

#2. Square-Well Potential
def square_well(r, epsilon=0.01, sigma=3.4, range=1.5):
    """Calculate the square well potential. 
    Parameters:
        r: distance between atoms. 
        epsilon: depth of potential well. 
        sigma: particle diameter
        range: range of the well. 
    Returns:
        float: Potential energy V(r)
    """
    if r < sigma:
        return 1000 # Approximate infinity as 1000
    elif sigma <= r < range * sigma:
        return -epsilon  # Well depth
    else:
        return 0  # Zero potential for distances beyond the well

#3. Lennard Jones Potential (from 2-1)
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    """ Calculates the Lennard Jones potential.
    Parameters:
        r: distance between the atoms.
        epsilon: depth of potential well. 
        sigma: distance at which potential is 0. 
    Returns:
        float: Potential energy V(r)
    """
    LJ = (4 * epsilon) * ((sigma / r)**12 - (sigma / r)**6)
    return LJ

#For integration
num_points = 1000
r_values = np.linspace(1e-3, 5 * sigma, num_points)

#Caclulating potentials
hard_sphere_potential = np.array([hard_sphere(r) for r in r_values])
square_well_potential = np.array([square_well(r) for r in r_values])
lennard_jones_potential = np.array([lennard_jones(r) for r in r_values])

# Compute B2V for each potential using trapezoidal integration at T=100 K
def compute_B2V(potential_function, temp):
    """Computes the second virial coefficient at constant volume. 
    Parameters:
        potential_function: Defines which type of potential equation will be used. 
    Returns:
        float: B2V(T)
    """
    denom = 1 / (k_b * temp)
    integrand = np.exp(-potential_function * denom) - 1
    B2V = (-2 * np.pi * N_a) * trapezoid(integrand * r_values**2, r_values)

    return B2V

B2V_hard_sphere = compute_B2V(hard_sphere_potential, temp=100)
print(B2V_hard_sphere)
B2V_square_well = compute_B2V(square_well_potential, temp =100)
print(B2V_square_well)
B2V_lennard_jones = compute_B2V(lennard_jones_potential, temp =100)
print(B2V_lennard_jones)

# Range of temperatures
temperatures = np.linspace(100, 800, 100)
B2V_hard_sphere = []
B2V_square_well = []
B2V_lennard_jones = []

# Compute B2V for each temperature in range
for temp in temperatures:
    B2V_hard_sphere.append(compute_B2V(hard_sphere_potential, temp))
    B2V_square_well.append(compute_B2V(square_well_potential, temp))
    B2V_lennard_jones.append(compute_B2V(lennard_jones_potential, temp))

# Convert to numpy arrays for plotting
B2V_hard_sphere = np.array(B2V_hard_sphere)
B2V_square_well = np.array(B2V_square_well)
B2V_lennard_jones = np.array(B2V_lennard_jones)

# Plotting
plt.plot(temperatures, B2V_hard_sphere, label='Hard Sphere')
plt.plot(temperatures, B2V_square_well, label='Square Well')
plt.plot(temperatures, B2V_lennard_jones, label='Lennard-Jones')
plt.axhline(0, color='black', label="$B_{2V} = 0$")
plt.title('Second Virial Coefficient $B_{2V}$ as a Function of Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('$B_{2V}$ ($Å^3/mol$)')
plt.legend()
plt.grid()
plt.xlim(100, 800)
plt.show()
plt.savefig("B2V_Temp")
