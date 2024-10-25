import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Constants
np.random.seed(42)
a_o = 1  # Bohr radius
R = 2  # Separation distance
L = 20  # Bounds for random sampling

def psi_2p_z(x, y, z):
    """
    Calculate the value of psi for the 2p_z orbital of hydrogen at a given point.
    Parameters:
    x: float, The x-coordinate of the point.
    y: float, The y-coordinate of the point.
    z: float, The z-coordinate of the point.
    Returns:
    psi: float, The value of the psi for 2p_z orbital at the given point.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    cos_theta = (z / r)
    
    # 2p_z orbital wave function
    psi = (1 / (4 * np.sqrt(2 * np.pi) * a_o**(3/2))) * ((r / a_o) * np.exp(-r / (2 * a_o))) * cos_theta

    return psi

# Set the integration limits
L = 20

# Loop over the number of points to sample
def random_sampling(n_points_list, L, R):
    '''Calculates the overlap integral using random sampling
    Parameters:
    n_points_list (nparray): array of points of sizes to be sampled over
    L (float): Integration limits
    R (float): Separation distance
    Returns:
    random_averages (nparray): array of overlap integral for every value in n_points)list
    '''
    random_averages = []
    for n_points in n_points_list:
        x = np.random.uniform(-L, L, n_points)
        y = np.random.uniform(-L, L, n_points)
        z = np.random.uniform(-L, L, n_points)
    
        # Calculate the integrand for the overlap integral
        integrand = psi_2p_z(x, y, (z + R/2)) * psi_2p_z(x, y, (z - R/2))  # Separation of 2 atomic units
    
        # Estimate the integral
        integral = (2*L)**3 * np.mean(integrand)
        random_averages.append(integral)
    return random_averages

x = 0
y = 0
z = np.linspace(-20, 20, 1000)
integrand = psi_2p_z(x, y, z+R/2) * psi_2p_z(x, y, z-R/2)
importance_sampling = norm.pdf(z)

def important_sampling(n_points_list, R):
    '''Calculates the overlap integral using normal distribution importance sampling
    Parameters:
    n_points_list (nparray): array of points of sizes to be sampled over
    R (float): Separation distance
    Returns:
    imp_averages (nparray): array of overlap integral for every value in n_points list
    '''
    imp_averages = []
    for n_points in n_points_list:
        x = norm.rvs(size=n_points, scale=3, loc=0)
        y = norm.rvs(size=n_points, scale=3, loc=0)
        z = norm.rvs(size=n_points, scale=3, loc=0)
        numer = psi_2p_z(x, y, (z+ (R/2))) * psi_2p_z(x, y, (z-(R/2)))
        denom = norm.pdf(x, scale=3, loc=0) * norm.pdf(y, scale=3, loc=0) * norm.pdf(z, scale=3, loc=0)
        integral = np.mean(numer/denom)
        imp_averages.append(integral)
    return imp_averages
