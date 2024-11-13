# Script for model 1 of the falling cat - Stable's model

### Importing libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### Parameters
R = 0.25 # radius of cat



### Coordinate system - spherical (with R fixed)



### Example of using Euler's method for solving ODEs

### Function: f(x) = dx/dt = ax
def f(t, x, a):
  """
  in:
  t: T dimensional ndarray such that ts[i]>0 and ts[i+1]>ts[i]
  x: TxN dimensional ndarray
  a: parameter

  out:
  a*x: TxN dimensional ndarray multiplied by parameter a
  """

  return np.array([a*x])

### Euler's method
def euler(f, x0, t):
  """
  in:
  f: python function which takes a N dimensional ndarray and returns an N dimensional ndarray
  x0: N dimensional ndarray
  t: T dimensional ndarray such that ts[i]>0 and ts[i+1]>ts[i]

  out:
  x: TxN dimensional ndarray
  """

  x = np.empty([len(t), len(x0)])
  x[0,:] = x0

  # loop over all time points in t array to update x
  for i in range(len(t)-1):

    h = t[i+1] - t[i] # the step size
    x[i+1,:] = x[i,:] + f(t[i], x[i,:], a) * h # update x using the Euler method

  return x

### Do something like this for the model :) TODO


### Plotting TODO
