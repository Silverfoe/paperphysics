# Author: Jacob Wyrozebski
# See provided link for details/documentation: https://docs.google.com/document/d/1LLWHzmtzaEDh65zR5MCt7fNyIwnT_0YYOere-0bIR3o/view?usp=sharing

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def compute_probability():
    print("\n--- Ultimate Reverse Order Paper Probability Calculator ---\n")

    try:
        # Basic Inputs (User Provides These)
        x = int(input("Enter the number of papers (must be ≥ 1): "))
        if x < 1:
            raise ValueError("Number of papers must be at least 1.")

        delta_rho = float(input("Enter air density fluctuation (∆ρ in kg/m³, e.g., 0.12): "))
        velocity = float(input("Enter falling velocity (m/s, e.g., 3.5): "))
        height = float(input("Enter drop height (m, e.g., 2.0): "))
        if velocity <= 0 or height <= 0:
            raise ValueError("Velocity and height must be greater than zero.")

        # Fundamental Constants (Fixed)
        Cd = 1.2  # Drag coefficient for flat paper
        Cl = 1.0  # Lift coefficient
        rho_air = 1.225  # Standard air density (kg/m³)
        A = 0.0625  # Paper surface area (m²)
        m = 0.0025  # Mass of a standard A4 paper (kg)
        g = 9.81  # Gravity (m/s²)
        t0 = 1.0  # Characteristic time scaling constant
        P_swap = 0.10  # 10% chance of swapping positions

        # Predefined Advanced Physics Factors
        quantum_fluctuations = 10**-25  # Quantum uncertainty
        dark_matter_interaction = 10**-22  # Dark matter gravity effects
        cosmic_entropy_disturbance = 10**-30  # Universe entropy disturbance
        multiversal_fluctuations = 10**-20  # Effects from parallel universes
        frame_dragging_effect = 10**-15  # Spacetime warping
        holographic_information_leak = 10**-18  # Holographic principle effects
        extra_dimensional_gravity = 10**-17  # Gravity from higher dimensions
        tachyonic_randomness = 10**-19  # Tachyonic influence
        hawking_radiation_entropy = 10**-35  # Quantum black hole effects
        kolmogorov_turbulence = 10**-10  # Computational turbulence complexity
        butterfly_effect_amplification = 10**-12  # Chaos theory sensitivity

        # Compute forces
        A_eff = A * 0.75  # Average effective area due to tilt
        F_d = 0.5 * Cd * rho_air * A_eff * velocity**2  # Drag force
        F_L = 0.5 * Cl * rho_air * A * velocity**2  # Lift force

        # Compute fall time
        v_term = math.sqrt((2 * m * g) / (rho_air * Cd * A))  # Terminal velocity
        t_fall = height / v_term  

        # Full Probability Disruption Function (Combining All Physics)
        D = (delta_rho + F_d + F_L +
             quantum_fluctuations + dark_matter_interaction + cosmic_entropy_disturbance +
             multiversal_fluctuations + frame_dragging_effect + holographic_information_leak +
             extra_dimensional_gravity + tachyonic_randomness + hawking_radiation_entropy +
             kolmogorov_turbulence + butterfly_effect_amplification)

        # Compute Final Probability
        try:
            probability = (math.exp(-D) * (1 - P_swap) ** (x - 1) * math.exp(-t_fall / t0)) / math.factorial(x)
        except OverflowError:
            probability = 0.0  # If factorial is too large, probability is effectively zero

        # Convert probability to scientific notation without loss of precision
        probability_percentage = probability * 100
        print(f"\nProbability of papers landing in reverse order: {probability_percentage:.50f} %\n")

        # Generate 3D Graphs
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
                    D = delta_rho + F_d + F_L + sum(10**-i for i in range(10, 40, 5))
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

# Run the optimized program
compute_probability()
