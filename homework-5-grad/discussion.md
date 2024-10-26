# Convergence for Diagonal and Off Diagonal KE Matrices

First diagonal matrices were examined. For random sampling, the value of Kii converges to ~0.50 at ~10^7 points sampled. On the other hand, Kii converges to the same value of 0.5 at ~10^5 points sampled, 2 orders of magnitude quicker. Both of these graphs are plotted on the same y scale to show how much more accurate importance sampling is in this case. 

Next off-diagonal matrices were examed, setting R equal to 1.4 and changing the z-axis. FOr radom sampling, Kij converges to ~2.1 similarly at 10^7 points sampled. Again for importance distrubution, the value converges to 2.1 in ~10^5 points sampled, 2 orders fo magnitude quicker. 

Again, my computer seems to run into an error with the kernal crashing when trying to calculate 10^8 points samples. To gather more data, 9^8 points sampled was collected. Additionally, when running importance sampling the kernel has problems when printing the integral value for 9^8 points sampled; however, it was okay to show in the graph. 