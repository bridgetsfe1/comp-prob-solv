# Difference in OLS and Optimization

In both cases, using OLS linear regression and minimizing least squares, the same results are returned for the values of a and b. This is as expected because the goals of both methods are to return the best line of fit, jsut using different methods. The only difference between the results is that a confidence interval was calucalted for the OLS, thus giving the value of slope and intercept a range for the OLS methods nto seen in the minimization of least squares. 

## Different Approches: Linear regression

In cases where the data is easily intepretable and only 2 variables are being tested, it is generally the case where one method would not be preferred over the other. Optimization however may take a shorter amount of time to code, but if all functions used in OLS are pre-made there would be little to no difference. 
When the data becomes more convoluted, optimization may be a better route. The minimization function can minimize a large amount of scalar values based on a function that is defined, which may take longer or be more confusing for an OLS model to do. 
In general, both lines of code are aiming to have the same result (minimizing least squares); however, the minimize function has a lot of pre-built ways to minimize the scalars which has to be done by hand for the OLS method. 