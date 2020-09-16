import matplotlib.pyplot as plt
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
# ---------------------- default figure setting ends---------------#


def ode_solver(expression, steps, step_size, initial_value, z0, tau=None):
    """This numerically the ODE using 4th order Runge-Kutta method
        ref: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#

    Parameters:
        expression (func): f(t,y) using the notation in ref above
        steps (int): Total number of steps needed to be integrated
        stepsize (float): The size of each step in t, i.e. Δ
        initial value (float): y(t_0)
        t0 (float): t_0 initial time of the integration
        tau (float): extra parameter for the function expression

    Returns:
        initial value (float): a generator for value of y at each step 
                                of integration
    """
    for i in range(steps):
        k1 = expression(z0, initial_value)
        k2 = expression(z0 + step_size/2, initial_value + step_size * k1 /2)
        k3 = expression(z0 + step_size/2, initial_value + step_size * k2 /2)
        k4 = expression(z0 + step_size, initial_value + step_size * k3)
        initial_value += step_size * (k1 + 2 * k2 + 2 * k3 + k4) /6
        z0 += step_size
        yield initial_value

# expression for the RHS of the differential equation
def dinten_dz(z, inten):
    return 5/(1 + inten/5) * inten

# 5 segments, define initial conditions
z1 = np.array([10** (-2) * i for i in range(5+1)])
I1 = np.zeros([5+1])
I1[0] = 10
i = 1
for tmp_res in ode_solver(dinten_dz, 5, 10** (-2), 10, 0, 1):
    I1[i] = tmp_res
    i += 1
print(I1[-1])

# 1000 segments, define initial conditions
z2 = np.array([0.05 * 10** (-3) * i for i in range(1000+1)])
I2 = np.zeros([1000+1])
I2[0] = 10
i = 1
for tmp_res in ode_solver(dinten_dz, 1000, 0.05 * 10** (-3), 10, 0):
    I2[i] = tmp_res
    i += 1
print(I2[-1])

# plot the results
fig, ax1 = plt.subplots()

plt.xlim(-0.005,0.055)
ax1.set_xlabel(r'$z/m$')
ax1.set_ylabel(r'$Intensity/(W/m^2)$') # , color=color)
ax1.plot(z2, I2, 'o',
            # linestyle='solid',
            ms=1, color='tab:red',
            label=r'$\textrm{1000 segments}$')
ax1.plot(z1, I1, 'o',
            # linestyle='solid',
            ms=5, color='tab:blue',
            label=r'$\textrm{5 segments}$')

ax1.legend() # loc='center right') 
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig("q4a.png")