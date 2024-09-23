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

initial_r = 4 

LJ_argon = minimize(lennard_jones, initial_r)

optimal_distance = LJ_argon.x[0]
print(f"Optimal distance: {optimal_distance:.2f} Å")

r_values = np.linspace(3, 6, 100)
V_values = lennard_jones(r_values)

# Plotting Potential Well and Optimization Value
plt.plot(r_values, V_values, label='Lennard-Jones Potential')
plt.axvline(optimal_distance, color='red',linestyle='--', label=f'Equilibrium Distance: {optimal_distance:.2f} Å')
plt.axhline(y=0, color='black')


# Labeling Graph
plt.title('Lennard-Jones Potential between 2 Argon Atoms')
plt.xlabel('Distance (Å)')
plt.ylabel('V(r)')
plt.xlim(3, 6)
plt.ylim(-.015, .03) 
plt.legend()
plt.grid()
plt.show()
plt.savefig("LJ_Argon_2.png")