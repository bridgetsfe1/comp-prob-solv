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
# Potential Matrix (Harmonic and Anharmonic)
def harmonic_matrix(grid):
    V_harmonic = 0.5 * m * omega**2 * (grid**2)
    return np.diag(V_harmonic) 

def anharmonic_matrix(grid):
    x_0 = 0  # You can change this if needed
    anh_exp = np.exp(-beta * (grid - x_0))
    V_anharmonic = D * (1 - anh_exp)**2
    return np.diag(V_anharmonic)  

#Laplacian Matrix
diagonal = -2 * np.ones(n) 
off_diagonal = np.ones(n-1) 

laplacian = (1/dx**2)*(np.diag(diagonal, k=0) + np.diag(off_diagonal, k=1) + np.diag(off_diagonal, k=-1))

#Hamiltonian Matrix
harmonic_V = harmonic_matrix(grid)
def hamil():
    hamiltonian_matrix = (-((h_bar**2) / (2 * m)) * laplacian) + harmonic_V
    return hamiltonian_matrix

anharmonic_V = anharmonic_matrix(grid)
def an_hamil():
    hamiltonian_matrix_morse = -((h_bar**2) / (2 * m)) * laplacian + anharmonic_V
    return hamiltonian_matrix_morse

#Compute eigenvalues and eigenvectors
eigenvalues_h, eigenvectors_h = np.linalg.eigh(hamil())
sort = np.argsort(eigenvalues_h)
eigenvalues_h = eigenvalues_h[sort]

eigenvalues_ah, eigenvectors_ah = np.linalg.eigh(an_hamil())
sort = np.argsort(eigenvalues_ah)
eigenvalues_ah = eigenvalues_ah[sort]

#Harmonic Graph                                          
for i in range(10):                               
    wavefunction_h= -eigenvectors_h[:, i]                                  
    wavefunction_h= wavefunction_h/(np.sqrt(np.sum(wavefunction_h**2)*dx))    #Normalize Wavefunction
    y=2*wavefunction_h+eigenvalues_h[i]                                     
    plt.plot(grid,y, label=f"n= $\psi_{i+1}$")                               
    plt.axvline(x=L/2)                  
    plt.axvline(x=-L/2)
    plt.hlines(eigenvalues_h[i],-L/2, L/2)               
plt.title("First Ten Wavefunctions with Energy Levels")
plt.xlabel("Length")
plt.ylabel("Energy (amu)")
plt.legend()
plt.savefig("First_10_harmonic")

#Anharmonic Graph
for i in range(10):                               
    wavefunction_ah= -eigenvectors_ah[:, i]                                  
    wavefunction_ah= wavefunction_ah/(np.sqrt(np.sum(wavefunction_ah**2)*dx))   
    y=2*wavefunction_ah+eigenvalues_ah[i]                                     
    plt.plot(grid,y, label=f"n= $\psi_{i+1}$")                               
    plt.axvline(x=L/2)                  
    plt.axvline(x=-L/2)
    plt.hlines(eigenvalues_ah[i],-L/2, L/2)               
plt.title("First Ten Wavefunctions with Energy Levels")
plt.xlabel("Length")
plt.ylabel("Energy (amu)")
plt.legend()
plt.savefig("First_10_anharmonic")

