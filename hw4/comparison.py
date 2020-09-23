import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib as plot
plot.use('Agg')
import numpy as np 
from scipy.integrate import odeint
# import helpfunction as hf

# ---------------------- default figure setting ---------------#
# plot.rc('text', usetex=True)
# plot.rc('font', family='cmunrm')
# plot.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plot.rcParams['mathtext.fontset'] = 'cm'
plot.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
plot.rcParams['xtick.major.size'] = 6
plot.rcParams['xtick.major.width'] = 1
plot.rcParams['xtick.minor.size'] = 3
plot.rcParams['xtick.minor.width'] = 1
plot.rcParams['ytick.major.size'] = 6
plot.rcParams['ytick.major.width'] = 1
plot.rcParams['ytick.minor.size'] = 3
plot.rcParams['ytick.minor.width'] = 1
plot.rcParams['axes.linewidth'] = 1.2
plot.rcParams['axes.edgecolor'] = 'black'
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['axes.facecolor'] = 'white'
plot.rcParams['text.color'] = 'black'
plot.rcParams['axes.labelcolor'] = 'black'
plot.rcParams['xtick.color'] = 'black'
plot.rcParams['ytick.color'] = 'black'

g1 = np.loadtxt('gain_euler.csv', delimiter=',')
g2 = np.loadtxt('gain_time_dependent_time=50ms.csv', delimiter=',')
g3 = np.loadtxt('gain_time_dependent_time=80ms.csv', delimiter=',')

fig, ax = plt.subplots()

xs1 = np.array([0.01 * i for i in range(1001)])
xs2 = np.array([0.01 * (i + 1) for i in range(1000)])
ax.plot(xs1, g1, label="steady state solution")
ax.plot(xs2, g2, label="time dependent solution after T=50ms")
ax.plot(xs2, g3, label="time dependent solution after T=80ms")
ax.legend()
fig.savefig("comparison.png")

ax.cla()
diff1 = np.zeros([1001])
diff2 = np.zeros([1001])
for i in range(1,1001):
    diff1[i] = g2[i-1] - g1[i]
    diff2[i] = g3[i-1] - g1[i]
ax.plot(xs1, diff1, label='50ms')
ax.plot(xs1, diff2, label='80ms')
ax.legend()
fig.savefig("difference.png")