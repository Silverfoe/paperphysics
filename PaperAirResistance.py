import math
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from mpl_toolkits.mplot3d import Axes3D # type: ignore

def compute_probability():
    print("\n--- Reverse Order Paper Probability Calculator ---\n")

    try:
        # User Inputs
        x = int(input("Enter the number of papers (must be ≥ 1): "))
        if x < 1:
            raise ValueError("Number of papers must be at least 1.")

        delta_rho = float(input("Enter air density fluctuation (∆ρ in kg/m³, e.g., 0.12): "))
        velocity = float(input("Enter falling velocity (m/s, e.g., 3.5): "))
        height = float(input("Enter drop height (m, e.g., 2.0): "))
        if velocity <= 0 or height <= 0:
            raise ValueError("Velocity and height must be greater than zero.")

        # Constants
        Cd = 1.2            # Drag coefficient
        Cl = 1.0            # Lift coefficient
        rho_air = 1.225     # Air density (kg/m³)
        A = 0.0625          # Paper surface area (m²)
        m = 0.0025          # Paper mass (kg)
        g = 9.81            # Gravity (m/s²)
        t0 = 1.0            # Time scaling constant
        P_swap = 0.10       # Probability of swapping positions

        # Scaling Constants (Estimates)
        k1 = 1.0  # Scaling factor for air density fluctuation
        k2 = 1.0  # Scaling factor for drag force
        k3 = 1.0  # Scaling factor for lift force

        # Compute Forces
        A_eff = A * 0.75                                    # Effective area
        F_d = 0.5 * Cd * rho_air * A_eff * velocity**2      # Drag force
        F_L = 0.5 * Cl * rho_air * A * velocity**2          # Lift force

        # Compute Fall Time
        v_term = math.sqrt((2 * m * g) / (rho_air * Cd * A))
        t_fall = height / v_term  

        # Compute Probability
        exponent = -(k1 * abs(delta_rho) + k2 * F_d + k3 * F_L)
        try:
            probability = (math.exp(exponent) * (1 - P_swap) ** (x - 1) * math.exp(-t_fall / t0)) / math.factorial(x)
        except OverflowError:
            probability = 0.0

        # Convert probability to scientific notation
        print(f"\nProbability of papers landing in reverse order: {probability:.3e} %\n")

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

    # Scaling Constants
    k1, k2, k3 = 1.0, 1.0, 1.0  # Use proper scaling factors

    for i, delta_rho in enumerate(delta_rho_values):
        for j, velocity in enumerate(velocity_values):
            for k, x in enumerate(x_values):
                try:
                    A_eff = 0.75 * 0.0625
                    F_d = 0.5 * 1.2 * 1.225 * A_eff * velocity**2
                    F_L = 0.5 * 1.0 * 1.225 * 0.0625 * velocity**2
                    v_term = math.sqrt((2 * 0.0025 * 9.81) / (1.225 * 1.2 * 0.0625))
                    t_fall = 2.0 / v_term

                    exponent = -(k1 * abs(delta_rho) + k2 * F_d + k3 * F_L)
                    probability = (math.exp(exponent) * (1 - 0.10) ** (x - 1) * math.exp(-t_fall / 1.0)) / math.factorial(x)
                except OverflowError:
                    probability = 0.0
                probabilities[i, j, k] = probability

    X, Y = np.meshgrid(x_values, delta_rho_values)
    Z = np.mean(probabilities, axis=1)

    # Convert Z to log10 only if probabilities are very small
    small_prob_threshold = 1e-25  # Adjust as needed
    if np.max(Z) < small_prob_threshold:
        Z = np.log10(Z + 1e-50)  # Apply log scale when values are very small
        z_label = "Log₁₀(Probability)"
    else:
        z_label = "Probability"

    # Plot Graph
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='plasma')
    ax.set_xlabel('Number of Papers (x)')
    ax.set_ylabel('Air Density Fluctuation (Δρ in kg/m³)')
    ax.set_zlabel(z_label)
    ax.set_title('Probability of Reverse Order Paper Landing')
    plt.show()

# Run Program
compute_probability()
