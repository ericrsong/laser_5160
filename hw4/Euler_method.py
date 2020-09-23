import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib as plot
plot.use('Agg')
import numpy as np 
# from scipy.integrate import odeint

def evolve(I_0, g, step_size_z=0.01):
    steps_z = g.shape[0]
    for j in range(steps_z):
        I_0 += step_size_z * (I_0 * g[j])
    return I_0

g1 = np.loadtxt('gain.csv', delimiter=',')[1:]
# print(g1.type())
g2 = np.loadtxt('gain_time_dependent_time=50ms.csv', delimiter=',')
# print(g2.type())
g3 = np.loadtxt('gain_time_dependent_time=80ms.csv', delimiter=',')

print("Steady state Euler Method yields")
print(evolve(10 ** (-3), g1))

print("50 ms Euler Method yields")
print(evolve(10 ** (-3), g2))

print("80 ms Euler Method yields")
print(evolve(10 ** (-3), g3))

