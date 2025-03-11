# Author: Jacob Wyrozebski
# See provided link for details/documentation: https://docs.google.com/document/d/1LLWHzmtzaEDh65zR5MCt7fNyIwnT_0YYOere-0bIR3o/view?usp=sharing

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def compute_probability():
    print("\n--- Reverse Order Paper Probability Calculator ---\n")
    
    try:
        x = int(input("Enter the number of papers (must be ≥ 1): "))
        if x < 1:
            raise ValueError("Number of papers must be at least 1.")

        delta_rho = float(input("Enter air density fluctuation (∆ρ in kg/m³, e.g., 0.12): "))
        velocity = float(input("Enter falling velocity (m/s, e.g., 3.5): "))
        height = float(input("Enter drop height (m, e.g., 2.0): "))
        if velocity <= 0 or height <= 0:
            raise ValueError("Velocity and height must be greater than zero.")

        # Constants
        Cd, Cl, rho_air, A, m, g, t0, P_swap = 1.2, 1.0, 1.225, 0.0625, 0.0025, 9.81, 1.0, 0.10
        A_eff = A * 0.75
        
        # Compute forces
        F_d = 0.5 * Cd * rho_air * A_eff * velocity**2
        F_L = 0.5 * Cl * rho_air * A * velocity**2
        
        # Compute fall time
        v_term = math.sqrt((2 * m * g) / (rho_air * Cd * A))
        t_fall = height / v_term
        
        # Disruption Function
        D = delta_rho + F_d + F_L
        
        # Compute Probability
        try:
            probability = (math.exp(-D) * (1 - P_swap) ** (x - 1) * math.exp(-t_fall / t0)) / math.factorial(x)
        except OverflowError:
            probability = 0.0

        print(f"\nProbability of papers landing in reverse order: {probability * 100:.10e} %\n")
        generate_3d_graphs(x, delta_rho, velocity)
    
    except ValueError as e:
        print(f"\nError: {e}. Please enter valid numerical values.\n")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}\n")

def generate_3d_graphs(x_value, delta_rho_value, velocity_value):
    x_values = np.arange(1, min(20, x_value+1))
    delta_rho_values = np.linspace(0, 0.5, 5)
    velocity_values = np.linspace(1, 10, 5)
    probabilities = np.zeros((len(delta_rho_values), len(velocity_values), len(x_values)))
    
    for i, delta_rho in enumerate(delta_rho_values):
        for j, velocity in enumerate(velocity_values):
            for k, x in enumerate(x_values):
                try:
                    A_eff = 0.75 * 0.0625
                    F_d = 0.5 * 1.2 * 1.225 * A_eff * velocity**2
                    F_L = 0.5 * 1.0 * 1.225 * 0.0625 * velocity**2
                    v_term = math.sqrt((2 * 0.0025 * 9.81) / (1.225 * 1.2 * 0.0625))
                    t_fall = 2.0 / v_term
                    D = delta_rho + F_d + F_L
                    probability = (math.exp(-D) * (1 - 0.10) ** (x - 1) * math.exp(-t_fall / 1.0)) / math.factorial(x)
                except OverflowError:
                    probability = 0.0
                probabilities[i, j, k] = probability
    
    X, Y = np.meshgrid(x_values, delta_rho_values)
    Z = np.log10(np.mean(probabilities, axis=1) + 1e-50)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='plasma')
    ax.set_xlabel('Number of Papers (x)')
    ax.set_ylabel('Air Density Fluctuation (Δρ in kg/m³)')
    ax.set_zlabel('Log10(Probability)')
    ax.set_title('Probability of Reverse Order Paper Landing')
    plt.show()

compute_probability()
