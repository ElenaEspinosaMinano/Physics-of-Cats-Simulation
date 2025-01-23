# Script for model 1 of the falling cat - Stable's model

### Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.ticker import MultipleLocator, FuncFormatter

### Default color cycle of matplotlib
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
fontsize=14
fontsize_title=18
fontsize_axes=14
fontsize_tick=fontsize
fontsizeSmall=14
fontsizeTiny=12

### Parameters
g = 9.81 # acceleration due to gravity (in m/s^2)
R = 0.1 # radius of cat (in meters)
m = 3 # mass of cat (in kg)
gamma = 2

# Define equations
# d(theta)/dt = theta_dot
# d(theta_dot)/dt = -g/R * sin(theta)

def f(t, y):
    theta = y[0] # theta
    theta_dot = y[1] # angular velocity

    dtheta_dt = theta_dot
    dtheta_dot_dt = g / R * np.sin(theta) - gamma * theta_dot
    
    return np.array([dtheta_dt, dtheta_dot_dt])

### Initial conditions
theta_0 = 0 # initial angle (in radians)
theta_dot_0 = 0.2 # initial angular velocity (in radians/s)

y0_list = [theta_0, theta_dot_0]
y0_array = np.array([theta_0, theta_dot_0])

t0 = 0.0
tf = 10.0
h = 0.001
N = (tf - t0) / h

t_array = np.linspace(t0, tf, int(N))  # time array
param = [g, R]  # parameters: g and R

### Solve using RK45 method
print(f"Solving DE with [theta, theta_dot] = [{theta_0} rad, {theta_dot_0} rad/s], [t0, tf] = [{t0}, {tf}], step size = {h}, no of steps = {N} ")
solution = solve_ivp(fun=f, t_span=[t0, tf], y0=[theta_0, theta_dot_0], method='RK45', t_eval=t_array)

# solution = euler(t_array, f, y0_array) # euler method doesn't work

### Get theta and theta_dot from the solution
theta = solution.y[0] # solution[:, 0]
theta_dot = solution.y[1]

# print(theta)
# print(theta_dot)

print(theta[-1])
print(theta_dot[-1])

z = R * np.cos(theta)

### Plot the results
fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)

# Height vs. time
axs[0].plot(t_array, z, label='Height (z)', color=colors[0])
axs[0].set_ylabel('Height, z\n (in meters)', fontsize=fontsize_axes)
axs[0].legend(fontsize=fontsizeTiny)
axs[0].grid(True, alpha=0.5)
axs[0].tick_params(axis='both', labelsize=fontsize_tick)

# Theta vs. time
axs[1].plot(t_array, theta, label='Angle (θ)', color=colors[1])
axs[1].set_ylabel('Angle, θ\n (in radians)', fontsize=fontsize_axes)
axs[1].legend(fontsize=fontsizeTiny)
axs[1].grid(True, alpha=0.5)
axs[1].tick_params(axis='both', labelsize=fontsize_tick)

### y-axis of theta plot as multiples of π
def pi_formatter(x, pos):
    return f"{x/np.pi:.0g}π" if x != 0 else "0"

axs[1].yaxis.set_major_locator(MultipleLocator(base=np.pi/2))
axs[1].yaxis.set_minor_locator(MultipleLocator(base=np.pi/4))
axs[1].yaxis.set_major_formatter(FuncFormatter(pi_formatter))

# Angular velocity vs. time
axs[2].plot(t_array, theta_dot, label='Angular Velocity (ω)', color=colors[2])
axs[2].set_xlabel('Time, t (in seconds)', fontsize=fontsize_axes)
axs[2].set_ylabel('Angular Velocity, ω\n (in radians/second)', fontsize=fontsize_axes)
axs[2].legend(fontsize=fontsizeTiny)
axs[2].grid(True, alpha=0.5)
axs[2].tick_params(axis='both', labelsize=fontsize_tick)

plt.tight_layout()
plt.savefig('../../Figures/'+'M1_via_angles_graphs.png', dpi=300, bbox_inches='tight')
plt.show()

### Plot the results
# plt.figure()
# # plt.plot(t_array, theta, label='theta (angle)')
# # plt.plot(t_array, theta_dot, label='theta_dot (angular velocity)')
# plt.plot(t_array, z, label='z (height)')
# plt.xlabel('Time, t (in seconds)')
# plt.ylabel('Height, z (in meters)')
# plt.legend()
# plt.show()


"""
    *** YAYYYYYY IT WORKS !!! ***
"""




### Could use Euler method to solve DE - doesn't work as DE seems to be stiff
def euler(t_array, f, y0):
    t_size = t_array.size
    y0_size = y0.size # number of equations (here: theta and theta_dot)

    Y = np.zeros((t_size, y0_size))
    
    Y[0,:] = y0 # initial values of theta and theta_dot

    # loop over time points and update Y using the Euler method
    for i in range(t_size - 1):
        dt = t_array[i+1] - t_array[i] # time step
        Y[i+1,:] = Y[i,:] + f(t_array[i], Y[i,:]) * dt

    return Y


"""
    INITIAL ROUGH CODE - doesn't work 
"""

"""import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### Parameters
g = 9.81 # acceleration due to gravity
R = 0.25 # radius of cat
m = 5 # mass of cat

# beta theta = first derivative of theta = d(theta)/dt = sin(theta) * g/R * t
def beta_theta(t, y, param):

    dtheta_dt = np.sin(y[0]) * g/R * t
    
    return np.array([dtheta_dt])

def beta_dot_theta(t, y, param):
    dbeta_dt = np.sin(y[0]) * g/R

    return np.array([dbeta_dt])

def euler(t_array, f, y0):
    t_size = t_array.size
    y0_size = y0.size # no of eqns

    Y = np.zeros((t_size, y0_size))
    
    Y[0,:] = y0 # the first values in the zero index of Y is set to array y0

    # loop over time points in t_array + update Y using Euler method
    for i in range(t_size - 1):
        
        dt = t_array[i+1] - t_array[i] # find the size of the time step  dt
        Y[i+1,:] = Y[i,:] + f(t_array[i], Y[i,:], param) * dt

    return Y

theta_ddot_0 = np.array([0.0, 0.0, 0.0]) # initial array of [theta_double_dot, phi_double_dot, z_double_dot]
#phi_ddot_0 = np.array([0.0, 0.0, 0.0]) # initial array of [theta_double_dot, phi_double_dot, z_double_dot]

t0 = 0.0 # initial value of t
tf = 100.0 # final value of t

h = 0.05 # step size
N = (tf-t0) / h

t_array = np.linspace(t0, tf, int(N+1)) # time array from t0 to tf with timestep h (note arange excludes endpoint, add h)
param = [g, R] # list of parameters: g and R

beta_new = euler(t_array, beta_dot_theta, theta_ddot_0) # obtain solution using appropiate model and initial conditions
theta_ddot = euler(t_array, beta_new, theta_ddot_0)

# Plot the results
plt.figure()
plt.plot(t_array, beta_new[:,0], label='theta_dot')
plt.plot(t_array, theta_ddot[:,0], label='theta_double_dot')
plt.xlabel('Time')
plt.ylabel('theta_dot')

plt.legend()
plt.show()"""
