# Overlap Integral Using Random Sampling versus Importance Sampling
The method of random sampling has shown the the overlap integral converges to approximately 0.73 at 10^5 points sampled. 

The method using importance sampling has shown a similar convergence to 0.73 at 10^4 points sampled, an order of magnitude quicker than random sampling. 

The method used to calculate importance sampling is normal distribution instead of exponential because the orbital shape fits a normal distribution better than a uniform distribution as used in random sampling or exponential distribution used in class for the 1s orbital. Choosing a sample that matches what is expected in the distribution thus leads to a reduced computational cost. 

# Overlap Integral versus Separation Distance
As the distance between R increases, the overlap integral decreases, which is obvious because the further theparticles are the less overlap they have. There is a change in sign where the overlap integral is negative because when the orbitals are too far from each other anti-bonding interactions over run bonding interactions because the two particles are out of phase. 