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

initial_r = 4                #initial guess for distance btwn atoms

LJ_argon = minimize(lennard_jones, initial_r)
optimal_distance = LJ_argon.x[0]

atom_dict_Ar = {'Ar1': np.array([0.0, 0.0, 0.0]),
    'Ar2': np.array([optimal_distance, 0.0, 0.0]),}

#Length fn from Homework 1-2def compute_bond_length(atom_dict, coord1, coord2):
def compute_bond_length(atom_dict, coord1, coord2):
    """ Calculates the distance in angstroms between two atoms using Cartesian coordinates.
    Parameters:
        dict: Dictionary where positions of atoms in a given molecule are stored. 
        coord1: atom key 1
        coord2: atom key 2
    Returns:
        float: the distance between the two points in angstroms
        Error if both atoms are not defined in the molecule's dictionary.
        Warning if the distance is longer than 2 Angstroms (longer than covalent bonds).
    """
    warning = " "
    if coord1 not in atom_dict or coord2 not in atom_dict:
        print("Error: Both coordinates must appear in the same molecule.")
    else:
        coordinate_1 = atom_dict[coord1]
        coordinate_2 = atom_dict[coord2]
        distance = math.sqrt((coordinate_2[0]-coordinate_1[0])**2 + (coordinate_2[1]-coordinate_1[1])**2 + (coordinate_2[2]-coordinate_1[2])**2)

        if distance > 2:
            warning = "Warning: The distance between the 2 atoms below is greater than 2 angstroms and is not a reasonable range for covalent bonds."
            return distance
        else:
            return distance

length_r = compute_bond_length(atom_dict_Ar, 'Ar1', 'Ar2')
print(length_r)

#Caclulate graph values
r_values = np.linspace(3, 6, 100)
V_values = lennard_jones(r_values)

# Plotting Potential Well, Optimization Value, and V(r)=0
plt.plot(r_values, V_values, label='Lennard-Jones Potential')
plt.axvline(optimal_distance, color='red',linestyle='--', label=f'Equilibrium Distance: {optimal_distance:.2f} Å')
plt.axhline(y=0, color='black')

# Labeling Graph
plt.title('Lennard-Jones Potential between 2 Argon Atoms')
plt.xlabel('Distance (Å)')
plt.ylabel('V(r) (arbitrary units)')
plt.xlim(3, 6)
plt.ylim(-.0125, .025) 
plt.legend()
plt.grid()
plt.show()
plt.savefig("LJ_Argon_2.png")

# Output to XYZ file
with open("argon_dimer.xyz", "w") as f:
    f.write("2\n")  # Number of atoms
    f.write("Optimal Argon Dimer\n")  #Molecule description
    for key in atom_dict_Ar:
        f.write(f"Ar  {atom_dict_Ar[key][0]:.3f}  {atom_dict_Ar[key][1]:.3f} {atom_dict_Ar[key][2]:.3f}\n")
