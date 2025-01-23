# Script for model 2 of the falling cat - Peano's model

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

k = 1 # some sort of spring constant between the cylinders

gamma = 2 # assume gamma_1 = gamma_2 = gamma = 2 -> damping coefficient??

gamma_1 = 0.2
gamma_2 = 0.2

tau_1 = 0.5 # internal torque of cylinder 1 (from cat!)
tau_2 = -0.5 # internal torque of cylinder 2 (from cat!)

I_1 = 1/2*m*R**2 # moment of inertia of cylinder 1
I_2 = 1/2*m*R**2 # moment of inertia of cylinder 2

# Define equations
# d(theta_1)/dt = v_1
# d(v_1)/dt = 1/I_1 * (tau_1 - gamma_1*v_1 - k*(theta_1 - theta_2))
# d(theta_2)/dt = v_2
# d(v_2)/dt = 1/I_2 * (tau_2 - gamma_2*v_2 - k*(theta_2 - theta_1))

def f(t, y):
    theta_1 = y[0]
    v_1 = y[1]
    theta_2 = y[2]
    v_2 = y[3]

    dtheta_1_dt = v_1
    dv_1_dt = 1/I_1 * (tau_1 - gamma_1*v_1 - k*(theta_1 - theta_2))
    dtheta_2_dt = v_2
    dv_2_dt = 1/I_2 * (tau_2 - gamma_2*v_2 - k*(theta_2 - theta_1))
    
    return np.array([dtheta_1_dt, dv_1_dt, dtheta_2_dt, dv_2_dt])

def Runge_Kutta_4(times, f, y0):
    no_of_times = times.size
    no_of_eqns = y0.size

    Y = np.zeros((no_of_times, no_of_eqns))

    Y[0,:] = y0 # the first values in the zero index of Y is set to array y0
    
    # loop over time points and update Y using RK4
    for i in range(no_of_times-1):
        dt = times[i+1] - times[i] # find the size of the time step dt
        
        # calculate the four stage variables
        delta1 = f(times[i],Y[i,:])* dt
        delta2 = f(times[i] + 0.5 * dt, Y[i,:] + 0.5 * delta1)* dt
        delta3 = f(times[i] + 0.5 * dt, Y[i,:] + 0.5 * delta2) * dt
        delta4 = f( times[i] + dt,Y[i,:] + delta3) * dt

        # update Y using RK4 method
        Y[i+1,:] = Y[i,:] + (delta1 + 2.0 * delta2 + 2.0 * delta3 + delta4)/6.
    
    return Y

### Initial conditions
theta_1_0 = 0 # initial angle (in radians)
v_1_0 = 0.2 # initial angular velocity (in radians/s)
theta_2_0 = 0 # initial angle (in radians)
v_2_0 = 0.2 # initial angular velocity (in radians/s)

# y0_list = [theta_1_0, v_1_0, theta_2_0, v_2_0]
# y0_array = np.array([theta_1_0, v_1_0, theta_2_0, v_2_0])

t0 = 0.0
tf = 10.0
h = 0.001
N = (tf - t0) / h

t_array = np.linspace(t0, tf, int(N))  # time array

### Solve using RK45 method
# print(f"Solving DE with [theta, theta_dot] = [{theta_0} rad, {theta_dot_0} rad/s], [t0, tf] = [{t0}, {tf}], step size = {h}, no of steps = {N} ")
# solution = solve_ivp(fun=f, t_span=[t0, tf], y0=[theta_1_0, v_1_0, theta_2_0, v_2_0], method='RK45', t_eval=t_array)

### Get theta and v from the solution - using scipy.solve_ivp
# theta_1 = solution.y[0]
# v_1 = solution.y[1]
# theta_2 = solution.y[2]
# v_2 = solution.y[3]

### Get theta and v - using the RK4 function (Davide likes to see where it comes from?)
Y = Runge_Kutta_4(t_array, f, np.array([theta_1_0, v_1_0, theta_2_0, v_2_0]))
print(Y.shape)
theta_1 = Y[:,0]
v_1 = Y[:,1]
theta_2 = Y[:,2]
v_2 = Y[:,3]

### Should now calculate z_1 and z_2 and then get the average (ie. the z_CoM)

### plot results
fig, axs = plt.subplots(1, 2, figsize=(16, 6))

### theta vs time
axs[0].plot(t_array, theta_1, label=r'$\theta_1$', alpha=0.7, color=colors[0])
axs[0].plot(t_array, theta_2, label=r'$\theta_2$', alpha=0.7, color=colors[1])
axs[0].set_xlabel('Time (s)', fontsize=fontsize_axes)
axs[0].set_xlim([0, 10])
axs[0].set_ylabel('Angle, $\\theta$\n(in radians)', fontsize=fontsize_axes)
axs[0].set_ylim([-np.pi/2, np.pi/2])
axs[0].legend(fontsize=fontsizeSmall)
axs[0].grid(True, alpha=0.5)
axs[0].tick_params(axis='both', labelsize=fontsize_tick)

### y-axis of theta plot as multiples of π
def pi_formatter(x, pos):
    return f"{x/np.pi:.0g}π" if x != 0 else "0"

axs[0].yaxis.set_major_locator(MultipleLocator(base=np.pi/2))
axs[0].yaxis.set_minor_locator(MultipleLocator(base=np.pi/4))
axs[0].yaxis.set_major_formatter(FuncFormatter(pi_formatter))

### angular velocity vs time
axs[1].plot(t_array, v_1, label=r'$v_1$', alpha=0.7, color=colors[2])
axs[1].plot(t_array, v_2, label=r'$v_2$', alpha=0.7, color=colors[3])
axs[1].set_xlabel('Time (s)', fontsize=fontsize_axes)
axs[1].set_ylabel('Angular Velocity, $\omega$\n(in radians/second)', fontsize=fontsize_axes)
axs[1].legend(fontsize=fontsizeSmall)
axs[1].grid(True, alpha=0.5)
axs[1].tick_params(axis='both', labelsize=fontsize_tick)

plt.tight_layout()
# plt.savefig('../../Figures/'+'Side_by_side_angles_and_velocities.png', dpi=300, bbox_inches='tight')
plt.show()


### alternative way to plot
plt.plot(t_array, theta_1, label=r'$\theta_1$', alpha=0.7)
plt.plot(t_array, theta_2, label=r'$\theta_2$', alpha=0.7)
plt.xlabel('Time (s)', fontsize=fontsize)
plt.ylabel('Angle (rad)', fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.grid(True, alpha=0.5)
plt.tick_params(axis='both', labelsize=fontsize_tick)
plt.show()

plt.plot(t_array, v_1, label=r'$v_1$', alpha=0.7)
plt.plot(t_array, v_2, label=r'$v_2$', alpha=0.7)
plt.xlabel('Time (s)', fontsize=fontsize)
plt.ylabel('Angular velocity (rad/s)', fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.grid(True, alpha=0.5)
plt.tick_params(axis='both', labelsize=fontsize_tick)
plt.show()