import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters specific to the cat's movement
m = 4.0  # Mass of the cat, approximately 4 kg
g = 9.81  # Gravitational acceleration (Earth's gravity)
k = 0.5  # Elastic constant of the potential field to simulate agility
R = 0.3  # Radius of equilibrium position (cat's typical movement range)
gamma = 0.05  # Damping coefficient to reduce oscillations, simulating quick reactions

# Define the gradient of the potential field (i.e., -∇V)
def grad_V(x, y, z):
    return np.array([
        k * x,        # Partial derivative of V with respect to x
        k * y,        # Partial derivative of V with respect to y
        m * g + k * z  # Partial derivative of V with respect to z
    ])

# Define the system of ODEs
def ode_system(t, state):
    x, y, z, vx, vy, vz = state
    fx, fy, fz = -grad_V(x, y, z) - gamma * np.array([vx, vy, vz])  # Force with damping term
    
    # Return the derivatives [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
    return [vx, vy, vz, fx/m, fy/m, fz/m]

# Initial conditions
x0, y0, z0 = 0.1, 0.1, 0.1  # Initial position
vx0, vy0, vz0 = 0, 0, 0     # Initial velocity
initial_state = [x0, y0, z0, vx0, vy0, vz0]

# Time span for the simulation
t_span = (0, 600)  # Time range for the simulation (short duration)
t_eval = np.linspace(t_span[0], t_span[1], 500)  # Time points for evaluation

# Solve the ODE
solution = solve_ivp(ode_system, t_span, initial_state, t_eval=t_eval)

# Extract the solutions
x, y, z = solution.y[0], solution.y[1], solution.y[2]

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(t_eval, z, label="z(t)")  # Plot the vertical position (z) over time
plt.xlabel("Time (t) [s]")
plt.ylabel("Position [m]")
plt.legend()
plt.title("Cat Position vs. Time")
plt.grid()
plt.show()

