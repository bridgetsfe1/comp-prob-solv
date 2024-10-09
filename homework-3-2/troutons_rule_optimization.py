import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import seaborn as sns

file_path = "/workspaces/comp-prob-solv/homework-3-1/trouton.csv"
df = pd.read_csv(file_path)

#Convert for kcal to joules
df['H_v (J/mol)'] = df['H_v (kcal/mol)'] * 4184

def objective(params, x, y):
    '''Objective function to calculate the sum of squared residuals for fitting a linear model.
    Parameters:
    - params: A list or array containing the parameters [a, b] to optimize.
    - x: An array of dependent variables.
    - y: An array of independent variables.

    Returns:
    float: Sum of squared residuals.
    '''
    a, b = params
    residuals = y - (a * x + b)
    return np.sum(residuals**2)

#Initial Guess
a_initial, b_initial = np.polyfit(df['T_B (K)'], df['H_v (J/mol)'], 1)
initial_params = [a_initial, b_initial]

# Perform the optimization
result = minimize(objective, initial_params, args=(df['T_B (K)'].values, df['H_v (J/mol)'].values))
print(result)

# Extract optimized parameters
a_opt, b_opt = result.x

# Create a new set of predicted H_v values for plotting
H_v_pred = a_opt * df['T_B (K)'] + b_opt


plt.figure(constrained_layout=True)
# Scatter plot with points color-coded by class
sns.scatterplot(x=df["T_B (K)"], y=df["H_v (J/mol)"], hue=df["Class"], palette='bright')

#Add optimized line of fit
plt.plot(df['T_B (K)'], H_v_pred, color='black', label='$H_v = a*T_b + b$')

#Adding a and b values on graph
plt.text(0.25, 0.60, f'a: {a_opt:.3f} J/mol*K\nb: {b_opt/1000:.3f} KJ',
         fontsize=10, ha='right', va='bottom', transform=plt.gca().transAxes)


# Labeling the plot
plt.title("Trouton's Rule Optimization")
plt.xlabel("Temperature $T_B$ (K)")
plt.ylabel("Enthalpy $H_v (J/mol)$")
plt.legend(loc = 'upper left')


plt.savefig("troutons_opt")