# Script for model 1 of the falling cat - Stable's model

### Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

### Parameters
g = 9.81 # acceleration due to gravity (in m/s^2)
R = 0.2 # radius of cat (in meters)
m = 3 # mass of cat (in kg)
gamma = 5

lam = 500000

### Define equations:
# Potential field
def potential_field(x, y, z):
    return m*g*z + lam*(x**2 + y**2 + z**2 - R**2)**2

# Gradient of potential field
def grad_potential_field(x, y, z):
    return np.array([
        4*lam*x*(x**2 + y**2 + z**2 - R**2),
        4*lam*y*(x**2 + y**2 + z**2 - R**2),
        m*g + 4*lam*z*(x**2 + y**2 + z**2 - R**2)
    ])

# Equation of motion
def cat_ode_mod(t, Y):
    x, y, z = Y[:3]
    vx, vy, vz = Y[3:]

    fx, fy, fz = ( -grad_potential_field(x, y, z) - gamma*np.array([vx, vy, vz]) ) / m

    return [vx, vy, vz, fx, fy, fz]

# Initial conditions - not needed (well can't use since y0 has to be 1D)
r0 = np.array([0.0, 0.0, float(R)])
v0 = np.array([0.0, 0.0, 0.0])

# Time span
t0 = 0.0
tf = 15.0
N = 1001 # Number of time points
T = np.linspace(t0, tf, N) # Time array

# Solve the ODE
solution = solve_ivp(fun=cat_ode_mod, t_span=[t0, tf], y0=[0.0, 0.0, float(R), 0.0, 0.0, 0.0], t_eval=T)

# Extract the solutions
x_t, y_t, z_t = solution.y[3], solution.y[4], solution.y[5]

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(T, x_t, label="x(t)")
# plt.plot(T, y_t, label="y(t)")
plt.plot(T, z_t, label="z(t)")

plt.xlabel("Time (t) [s]")
plt.ylabel("Position [m]")
plt.legend()
plt.title("Cat Position vs. Time")
plt.grid()
plt.show()

exit()





























# Define equations

def Runge_Kutta_4(times,f,y0):
    no_of_times = times.size                          # find out how many different times required
    no_of_eqns = y0.size                              # find out how many equations
    Y = np.zeros((no_of_times, no_of_eqns))           # an array of zeros for storing the Y values for any number of eqns at 
                                                      # at each required time point
    
    Y[0,:] = y0 # the first values in the zero index of Y is set to array y0

    for i in range(no_of_times-1):                      # loop over time points
        dt = times[i+1] - times[i]                      # find the size of the time step  dt
        delta1 = f(times[i],Y[i,:])* dt             # calculate the first stage variable
        delta2 = f(times[i] + 0.5 * dt, Y[i,:] + 0.5  * delta1)* dt # calculate the second stage variable
        delta3 = f(times[i] + 0.5 * dt, Y[i,:] + 0.5 *  delta2) * dt # calculate the third stage variable
        delta4 = f( times[i] + dt,Y[i,:] +  delta3) * dt              # calculate the fourth stage variable
        Y[i+1,:] = Y[i,:] +  (delta1 + 2.0 * delta2 + 2.0 * delta3 + delta4)/6.0  #calculate Y[i+1,:]
    return Y                                                               # return Y

# def repressilator(t,y,par):                  # the repressilator function with arguments (t, y, par)
#     alpha0 = par[0]                           # alpha0 is the first element of the par array
#     n = par[1]                                 # n is the second element of the the par array
#     beta = par[2]                              # beta is the third element of the par array
#     alpha = par[3]                             # alpha is the fourth element of the par array
#     # y consists of an array [m1,p1,m2,p2,m3,p3]          
#     dm1_dt = -y[0]  +alpha/(1.0 + y[5]**n) + alpha0        # calculate dm1/dt
#     dp1_dt = -beta * (y[1] - y[0])                         # calculate dp1/dt
#     dm2_dt = -y[2]+ alpha/(1 + y[1]**n) + alpha0           # calculate dm2/dt
#     dp2_dt= -beta * (y[3] - y[2])                          # calculate dp2/dt
#     dm3_dt = -y[4] + alpha/(1 + y[3]**n) + alpha0          # calculate dm3/dt
#     dp3_dt = -beta * (y[5] - y[4])                         # calculate dp3/dt
    
#     # return the array[dm1/dt,dp1/dt, dm2/dt,dp2/dt, dm3/dt, dp3/dt]
#     return np.array([dm1_dt,dp1_dt,dm2_dt,dp2_dt,dm3_dt,dp3_dt]) 

# y0 = np.array([0.0 , 1.0 , 0.0, 2.0, 0.0, 3.0])  # the initial array of values of [m1,p1,m2,p2,m3,p3]
# t0 = 0.0     # the initial value of t
# tf = 15.0    # the final value of t

# N2 = 1001    # set the number of points N2 to 1001
                    
# T = np.linspace(t0, tf, N2)             # a T array to hold N2 values from t0 to tf
# par =[1.0,1.0,5.0, 1000.0]          # an array of parameters
# Y4 = Runge_Kutta_4(T, repressilator,y0)     # use the Runge-Kutta 4 stage method for the repressilator with T and y0

"""
    Now aim to modify the above to solve the cat problem
"""

def potential_field(x, y, z):
    return m*g*z + lam*(x**2 + y**2 + z**2 - R**2)**2

def grad_potential_field(x, y, z):
    return np.array([
        4*lam*x*(x**2 + y**2 + z**2 - R**2),
        4*lam*y*(x**2 + y**2 + z**2 - R**2),
        m*g + 4*lam*z*(x**2 + y**2 + z**2 - R**2)
    ])

def velocity(x, y, z):
    vx = 0
    vy = 0
    vz = 0

    return np.array([vx, vy, vz])

def cat_ode(t, r):

    fx, fy, fz = ( -grad_potential_field(r[0], r[1], r[2]) - gamma*velocity(r[0], r[1], r[2]) ) / m

    return [fx, fy, fz]

r0 = np.array([0.0, 0.0, 0.0])
t0 = 0.0
tf = 15.0

N = 1001

T = np.linspace(t0, tf, N)

Y = Runge_Kutta_4(T, cat_ode, r0)

x = Y[:,0]
y = Y[:,1]
z = Y[:,2]

plt.figure(figsize=(8, 6))
plt.plot(T, x, label="x(t)")
plt.plot(T, y, label="y(t)")
plt.plot(T, z, label="z(t)")

plt.xlabel("Time (t) [s]")
plt.ylabel("Position [m]")
plt.legend()
plt.title("Cat Position vs. Time")
plt.grid()
plt.show()





"""

# Define the gradient of the potential field (i.e., -∇V)
def grad_V(x, y, z):
    return np.array([
                        2*lam*x,
                        2*lam*y,
                        m * g + 2*lam*z
                    ])

# Define the system of ODEs
def ode_system(t, state):
    x, y, z, vx, vy, vz = state
    fx, fy, fz = 1/m * (-grad_V(x, y, z) - gamma * np.array([vx, vy, vz]))  # Force with damping term
    
    # Return the derivatives [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
    return [vx, vy, vz, fx/m, fy/m, fz/m]

# Initial conditions
x0, y0, z0 = 0.2, 0.2, 0.2  # Initial position
vx0, vy0, vz0 = 0, 0, 0     # Initial velocity
initial_state = [x0, y0, z0, vx0, vy0, vz0]

# Time span for the simulation
t_span = (0, 200)  # Time range for the simulation (short duration)
t_eval = np.linspace(t_span[0], t_span[1], 100)  # Time points for evaluation

# Solve the ODE
solution = solve_ivp(ode_system, t_span, initial_state, t_eval=t_eval)

# Extract the solutions
x, y, z = solution.y[0], solution.y[1], solution.y[2]
# print(solution.y[1], solution.y[2])

"""

# # Plot the results
# plt.figure(figsize=(8, 6))
# plt.plot(t_eval, x, label="x(t)")  # Plot the vertical position (z) over time
# plt.plot(t_eval, y, label="y(t)")  # Plot the vertical position (z) over time
# # plt.plot(t_eval, z, label="z(t)")  # Plot the vertical position (z) over time

# plt.xlabel("Time (t) [s]")
# plt.ylabel("Position [m]")
# plt.legend()
# plt.title("Cat Position vs. Time")
# plt.grid()
# plt.show()

