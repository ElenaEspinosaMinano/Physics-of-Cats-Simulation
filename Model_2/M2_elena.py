# Script for model 2 of the falling cat - Peano's model

### Importing libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

### Parameters
g = 9.81 # acceleration due to gravity (in m/s^2)
R = 0.25 # radius of cat (in meters)
m = 5 # mass of cat (in kg)

