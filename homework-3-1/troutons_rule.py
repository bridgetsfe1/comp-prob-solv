import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t
import seaborn as sns

df = pd.read_csv('trouton.csv')

#Convert for kcal to joules
df['H_v (J/mol)'] = df['H_v (kcal/mol)'] * 4184

y = df["H_v (J/mol)"]
x = df["T_B (K)"]

#All functions below copied from in class excersises
def ols_slope(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean) ** 2)
    return numerator / denominator

def ols_intercept(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    slope = ols_slope(x, y)
    return y_mean - slope * x_mean

def ols(x, y):
    slope = ols_slope(x, y)
    intercept = ols_intercept(x, y)
    return slope, intercept

# Perform OLS regression
slope, intercept = ols(x, y)

# Calculate residuals
residuals = y - (intercept + slope * x)

# Calculate the sum of the squared residuals
def sse(residuals):
    return np.sum(residuals ** 2)

# Calculate the variance of the residuals
def variance(residuals):
    return sse(residuals) / (len(residuals) - 2)

# Calculate the standard error of the slope
def se_slope(x, residuals):
    numerator = variance(residuals)
    x_mean = np.mean(x)
    denominator = np.sum((x - x_mean) ** 2)
    return np.sqrt(numerator / denominator)

# Calculate the standard error of the intercept
def se_intercept(x, residuals):
    numerator = variance(residuals)
    x_mean = np.mean(x)
    denominator = len(x) * np.sum((x - x_mean) ** 2)
    return np.sqrt(numerator / denominator)

#Calculates standard error values
se_slope_val = se_slope(x, residuals)
se_intercept_val = se_intercept(x, residuals)


#Calculate confiedence intervals for slope and intercept
def confidence_interval_slope(x, residuals, confidence_level):
    # Calculate the standard error of the slope
    se = se_slope(x, residuals)

    # Calculate the critical t-value
    n_data_points = len(x)
    df = n_data_points - 2  # degrees of freedom
    alpha = 1 - confidence_level
    critical_t_value = t.ppf(1 - alpha/2, df)

    # Calculate the confidence interval
    return critical_t_value * se
print(f"slope: {slope:.3f} +/- {confidence_interval_slope(x, residuals, 0.95):.3f}")
CI_slope = confidence_interval_slope(x, residuals, 0.95)


def confidence_interval_intercept(x, residuals, confidence_level):
    # Calculate the standard error of the intercept
    se = se_intercept(x, residuals)

    # Calculate the critical t-value
    n_data_points = len(x)
    df = n_data_points - 2  # degrees of freedom
    alpha = 1 - confidence_level
    critical_t_value = t.ppf(1 - alpha/2, df)

    # Calculate the confidence interval
    return critical_t_value * se
print(f"intercept: {intercept:.3f} +/- {confidence_interval_intercept(x, residuals, 0.95):.3f}")
CI_int = confidence_interval_intercept(x, residuals, 0.95)


#Plotting results
plt.figure(constrained_layout=True)
x_range = np.linspace(x.min(), x.max(), 100)
y_range = intercept + slope * x_range
plt.plot(x_range, y_range, color='black', linewidth=2, label='$H_v = a*T_b + b$')

plt.fill_between(x, (slope - CI_slope) * x + (intercept - CI_int),
    (slope + CI_slope) * x + (intercept + CI_int),
    color='black', alpha=0.3, label='95% CI')

# Scatter plot with points color-coded by class
sns.scatterplot(x=x, y=y, hue=df["Class"], palette='bright')

plt.text(0.95, 0.05, f'a: {slope:.3f} ± {CI_slope:.3f} J/mol*K\nb: {intercept/1000:.3f} ± {CI_int/1000:.3f} KJ/mol',
         fontsize=10, ha='right', va='bottom', transform=plt.gca().transAxes)

plt.title("Trouton's Rule")
plt.xlabel('Temperature $T_b$ (K)')
plt.ylabel('Enthalpy $H_v$ (J/mol)')
plt.legend()
plt.grid()

plt.savefig("troutons_rule")
