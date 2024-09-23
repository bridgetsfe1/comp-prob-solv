import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize


def lennard_jones(r, epsilon = 0.01, sigma = 3.4):
    """ Calculates the Lennard Jones potential between 2 Argon atoms. 
    Parameters:
        r: distance between the atoms.
        epsilon: depth of potential well. (0.01 for 2 Argon atoms)
        sigma: distance at which potential is 0. (3.4 for 2 Argon atoms)
    Return:
        float: Potential energy V(r)
    """

    LJ = (4 * epsilon) * ((sigma / r)**12 - (sigma / r)**6)
    return LJ

def total_potential(positions, epsilon=0.01, sigma=3.4):
    """Calculates the total Lennard Jones potential if the Argon trimer. Let Ar_1 be at the origin Ar_2 be on the x-axis and Ar_3 be unknown). 
    Parameters:
        variables: r_12 (distance between Ar_1 and Ar_2) x_3 and y_3 (position of Ar_3)
        epsilon: depth of potential well. 
        sigma: distance at which potential is 0. 

    Return: 
        float: total potential energy between 3 atoms V(r) = V_12(r) + V_13(r) + V_23(r)
    """

    r_12, x_3, y_3 = positions                       #distance between Ar_1 and Ar_2 and x, y coordinates of Ar_3

     # Calculate distances between the atoms (trimer)
    r_13 = np.sqrt(x_3**2 + y_3**2)
    r_23 = np.sqrt((x_3 - r_12)**2 + y_3**2) 

    # Calculate potential energy for each pair
    V13 = lennard_jones(r_13, epsilon, sigma)
    V23 = lennard_jones(r_23, epsilon, sigma)
    V12 = lennard_jones(r_12, epsilon, sigma)   #Distance between Ar_1 and Ar_2 optimized to be 3.82 in previous fn
    
    return V13 + V23 + V12

#Optimization initial guess and minimization
initial_positions = ([4, 2, 2])
result = minimize(total_potential, initial_positions)

optimal_r12, optimal_x3, optimal_y3 = result.x
print(f"Optimal distance between atoms 1 and 2: {optimal_r12:.2f} Å")
#print(f"Optimal position for atom 3: ({optimal_x3:.2f}, {optimal_y3:.2f})")

# Calculate r13 and r23
optimal_r13 = np.sqrt(optimal_x3**2 + optimal_y3**2)  # Distance between atom 1 and atom 3
optimal_r23 = np.sqrt((optimal_x3 - optimal_r12)**2 + optimal_y3**2)  # Distance between atom 2 and atom 3
    
print(f"Optimal distance between atoms 1 and 3: {optimal_r13:.2f} Å")
print(f"Optimal distance between atoms 2 and 3: {optimal_r23:.2f} Å")
    
# Calculate angles between argon atoms
theta_1 = np.arccos((optimal_r12**2 + optimal_r13**2 - optimal_r23**2) / (2 * optimal_r12 * optimal_r13))
theta_2 = np.arccos((optimal_r12**2 + optimal_r23**2 - optimal_r13**2) / (2 * optimal_r12 * optimal_r23))
theta_3 = np.arccos((optimal_r13**2 + optimal_r23**2 - optimal_r12**2) / (2 * optimal_r13 * optimal_r23))

print(f"Optimal angle at atom 1 (degrees): {np.degrees(theta_1):.2f}°")
print(f"Optimal angle at atom 2 (degrees): {np.degrees(theta_2):.2f}°")
print(f"Optimal angle at atom 3 (degrees): {np.degrees(theta_3):.2f}°")

print("Based on the output of equal bond distances and angles the atoms are positioned to create an equilateral triangle.")

