import math
import numpy as np
import matplotlib.pyplot as plt

#Constants
h_bar = 1
m = 1
omega = 1
D = 10
beta = np.sqrt(1/(2*D))
L = 40
n = 2000
dx = (L)/(n-1)

#Discretizing the Space
grid = np.linspace(-L / 2, L / 2, n)

#Potential Matrix (Harmonic and Anharmoic)
def harmonic_matrix(grid):
    """Computes the harmonic potential
    Parameters:
        grid: the size of the space in which potentials will be calculated
    Returns:
        list: Harmonic potentials over the grid space. 
    """
    V_harmonic = 0.5 * m * omega**2 * (grid**2)
    return (V_harmonic) 

def anharmonic_matrix(grid):
    """Computes the anharmonic (Morse) potential
    Parameters:
        grid: the size of the space in which potentials will be calculated
    Returns:
        list: Harmonic potentials over the grid space. 
    """
    x_0 = 0  # You can change this if needed
    anh_exp = np.exp(-beta * (grid - x_0))
    V_anharmonic = D * (1 - anh_exp)**2
    return (V_anharmonic)  

#Laplacian Matrix
diagonal = -2 * np.ones(n) 
off_diagonal = np.ones(n-1) 

laplacian = (1/dx**2)*(np.diag(diagonal, k=0) + np.diag(off_diagonal, k=1) + np.diag(off_diagonal, k=-1))

#Hamiltonian Matrix
harmonic_V = np.diag(harmonic_matrix(grid))
def hamil():
    """Computes the harmonic Hamiltonian
    """
    hamiltonian_matrix = (-((h_bar**2) / (2 * m)) * laplacian) + harmonic_V
    return hamiltonian_matrix

anharmonic_V = np.diag(anharmonic_matrix(grid))
def an_hamil():
    """Computes the anharmonic Hamiltonian
    """
    hamiltonian_matrix_morse = -((h_bar**2) / (2 * m)) * laplacian + anharmonic_V
    return hamiltonian_matrix_morse

#Compute and sort eigenvalues and eigenvectors (energies)
eigenvalues_h, eigenvectors_h = np.linalg.eigh(hamil())
sort = np.argsort(eigenvalues_h)
eigenvalues_h = eigenvalues_h[sort]

eigenvalues_ah, eigenvectors_ah = np.linalg.eigh(an_hamil())
sort = np.argsort(eigenvalues_ah)
eigenvalues_ah = eigenvalues_ah[sort]

#Harmonic Graph                                          
for i in range(10):                               
    wavefunction_h= eigenvectors_h[:, i]                                  
    y=2*wavefunction_h+eigenvalues_h[i]                                     
    plt.plot(grid,y, label=f"$\u03A8_{i+1}$: E = {eigenvalues_h[i]:.2f}")                             
    plt.axvline(x=L/2, color = 'black')                  
    plt.axvline(x=-L/2, color = 'black') 
plt.plot(grid, harmonic_matrix(grid), color = "grey", linestyle= "--", label = "Harmonic Potential")            
plt.title("First Ten Wavefunctions with Energy Levels")
plt.xlabel("Length (a.u.)")
plt.ylabel("Energy (amu)")
plt.ylim(0, 10)
plt.legend()
plt.savefig("First_10_harmonic")

#Anharmonic Graph
plt.clf()
for i in range(10):                               
    wavefunction_ah= eigenvectors_ah[:, i]                                  
    y=2*wavefunction_ah+eigenvalues_ah[i]                                     
    plt.plot(grid,y, label=f"$\u03A8_{i+1}$: E = {eigenvalues_ah[i]:.2f}")                               
    plt.axvline(x=L/2, color = 'black')                  
    plt.axvline(x=-L/2, color = 'black')   
plt.plot(grid, anharmonic_matrix(grid), color = "grey", linestyle= "--", label = "Morse Potential")                       
plt.title("First Ten Wavefunctions with Energy Levels")
plt.xlabel("Length (a.u.)")
plt.ylabel("Energy (amu)")
plt.ylim(0, 8)
plt.legend()
plt.savefig("First_10_anharmonic")

