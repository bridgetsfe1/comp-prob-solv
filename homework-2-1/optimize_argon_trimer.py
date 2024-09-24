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
        postions: r_12 (distance between Ar_1 and Ar_2) x_3 and y_3 (position of Ar_3)
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

atom_dict_Ar = {'Ar1': np.array([0.0, 0.0, 0.0]),
    'Ar2': np.array([optimal_r12, 0.0, 0.0]),
    'Ar3': np.array([optimal_x3, optimal_y3, 0.0])}

#Length fn from Homework 1-2
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

# Calculate bond lengths between argon atoms
#optimal_r12 is optimized in minimization fn
optimal_r13 = compute_bond_length(atom_dict_Ar, 'Ar1', 'Ar3')  # Distance between atom 1 and atom 3
optimal_r23 = compute_bond_length(atom_dict_Ar, 'Ar2', 'Ar3')  # Distance between atom 2 and atom 3
    
print(f"Optimal distance between atoms 1 and 2: {optimal_r12:.2f} Å")
print(f"Optimal distance between atoms 1 and 3: {optimal_r13:.2f} Å")
print(f"Optimal distance between atoms 2 and 3: {optimal_r23:.2f} Å")


#Angle fn from Homework 1-2
def compute_bond_angle(atom_dict, coord1, coord2, coord3):
    """ Calculate the angle between 3 atoms in a molecule in degrees using Cartesian coordinates. 
    Parameters:
        dict: Dictionary where positions of atoms in a given molecule are stored. 
        coord1: atom key 1
        coord2: atom key 2
        coord3: atom key 3
    Returns:
        Whether the bond angle is acute, right, or obtuse. 
        float: the bond angle between 3 atoms in degrees. 
        Error if all 3 atoms do not appear in the same dictionary.
    """
    if coord1 not in atom_dict or coord2 not in atom_dict or coord3 not in atom_dict:
        print("Error: All coordinates must appear in the same molecule.")
    else:
        coordinate_A = atom_dict[coord1]
        coordinate_B = atom_dict[coord2]
        coordinate_C = atom_dict[coord3]
        vector_BA = np.array([(coordinate_A[0]-coordinate_B[0]), (coordinate_A[1]-coordinate_B[1]), (coordinate_A[2]-coordinate_B[2])])
        vector_BC = np.array([(coordinate_C[0]-coordinate_B[0]), (coordinate_C[1]-coordinate_B[1]), (coordinate_C[2]-coordinate_B[2])])
        
        mag_AB = math.sqrt((vector_BA[0])**2 + (vector_BA[1])**2 + (vector_BA[2])**2)
        mag_BC = math.sqrt((vector_BC[0])**2 + (vector_BC[1])**2 + (vector_BC[2])**2)
        
        cos_angle = ((np.dot(vector_BA, vector_BC)) / (mag_AB * mag_BC))
        theta_rad = np.arccos(cos_angle)
        theta_deg = math.degrees(theta_rad)

        # if theta_deg == 90.00:
        #     print("The below bond angle is right.")
        # elif theta_deg > 90.00:
        #     print("The below bond angle is obtuse.")
        # else:
        #     print("The below bond angle is acute.")
    
    return theta_deg

# Calculate angles between argon atoms
theta_1 = compute_bond_angle(atom_dict_Ar, 'Ar2', 'Ar1', 'Ar3')
theta_2 = compute_bond_angle(atom_dict_Ar, 'Ar1', 'Ar2', 'Ar3')
theta_3 = compute_bond_angle(atom_dict_Ar, 'Ar1', 'Ar3', 'Ar2')

print(f"Optimal angle at atom 1 (degrees): {(theta_1):.2f}°")
print(f"Optimal angle at atom 2 (degrees): {(theta_2):.2f}°")
print(f"Optimal angle at atom 3 (degrees): {(theta_3):.2f}°")

print("Based on the output of equal bond distances and angles the atoms are positioned to create an equilateral triangle.")


# Output to XYZ file
with open("argon_trimer.xyz", "w") as f:
    f.write("3\n")  # Number of atoms
    f.write("Optimal Argon Trimer\n")  #Molecule description
    for key in atom_dict_Ar:
        f.write(f"Ar  {atom_dict_Ar[key][0]:.3f}  {atom_dict_Ar[key][1]:.3f}  {atom_dict_Ar[key][2]:.3f}\n")

