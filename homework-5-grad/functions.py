import numpy as np
from scipy.constants import pi
from scipy.stats import norm

def psi_1s(x, y, z, Z=1, a_o=1):
    """ Computes the 1s orbital wavefunction at a point (x,y,z)
    Parameters:
     x (float/nparray) : x-coordinate
     y (float/nparray) : y-coordinate 
     z (float/nparray) : z-coordinate 
     a_o (float): Bohr radius
     Returns:
     psi (float/nparray): Value of the 1s wavefunction at x,y,z
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    psi = (1/np.sqrt(pi * a_o**3)) * np.exp(-r / a_o)
    return psi

def laplacian_psi_1s(x, y, z, Z=1, a_o=1):
    """Computes the Laplacian of the hydrogen 1s orbital in a coordinate plane
    Parameters:
    x (float/nparray): x-coordinate
    y (float/nparray): y-coordinate
    z (float/nparray): z-coordinate
    Z (int): Atomic number
    a_o (float): Bohr radius
    Returns:
    float or array: The value of the Laplacian of the 1s orbital at the given point.
    """

    r = np.sqrt(x**2 + y**2 + z**2)
    psi = psi_1s(x, y, z, Z, a_o)

    d_psi_dr = -psi/ a_o         #First deriv
    d2_psi_dr2 = psi/ (a_o**2)   #Second deriv

    # Laplacian of the 1s orbital
    laplacian = d2_psi_dr2 + (2 / r) * d_psi_dr
    return laplacian

def random(L,N,R):
    """Computes the Monte Carlo integration using random sampling for the calculation of kinetic energy matrix of two 1s orbitals
    Parameters:
    L (int): Length of box
    N (int): Number of points for grid
    R (int): Separation distance between orbitals (0 being no separation)

    Returns:
    Kii (float): integration using random sampling
    """
    np.random.seed(42)

    #Constants
    Z = 1  
    a_o = 1  
    V = (2 * L)**3  
   
    x = np.random.uniform(-L, L, N)
    y = np.random.uniform(-L, L, N)
    z = np.random.uniform(-L, L, N)
    
    # Compute the integrand at each point
    psi = psi_1s(x, y, z+(R/2), Z, a_o)
    laplace = laplacian_psi_1s(x, y, z-(R/2), Z, a_o)
    integrand = -0.5 * psi * laplace
   
    # Estimate K_ii
    Kii = V * np.mean(integrand)                      

    return Kii

def importance(N,R):
    """Computes the Monte Carlo integration using normal distribution important sampling for the calculation of kinetic energy matrix of two 1s orbitals
    Parameters:
    N (int): Number of points
    R (int): Separation distance between the 2 1s orbitals
    Returns:
    Kii (float): integration using importance sampling
    """
    np.random.seed(42)

    #Constants
    Z = 1  # Atomic number for hydrogen
    a_o = 1  # Bohr radius in atomic units

    # Generate random points in the cubic region
    x = norm.rvs(size=N, scale=1)
    y = norm.rvs(size=N, scale=1)
    z = norm.rvs(size=N, scale=1)

    # Compute the integrand at each point
    psi = psi_1s(x, y, z+(R/2), Z, a_o)
    laplacian_psi = laplacian_psi_1s(x, y, z-(R/2), Z, a_o)
    numer = -0.5 * psi * laplacian_psi
    denom = norm.pdf(x) * norm.pdf(y) * norm.pdf(z)
    integrand=numer/denom

    # Monte Carlo estimation of the integral
    Kii = np.mean(integrand)            

    return Kii
 


