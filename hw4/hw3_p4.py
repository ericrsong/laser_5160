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
        average = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        gain = 5/(1 + initial_value/5)
        initial_value += step_size * average
        z0 += step_size
        yield (initial_value, gain)

# expression for the RHS of the differential equation
def dinten_dz(z, inten):
    return 5/(1 + inten/5) * inten


def solve(segments, total_length=10):
    step_size = total_length/segments
    # 1000 segments, define initial conditions
    # z2 = np.array([step_size * i for i in range(segments+1)])
    I2 = np.zeros([segments+1])
    I2[0] = 10 ** (-3)
    i = 1
    for (tmp_res, _) in ode_solver(dinten_dz, segments, step_size, I2[0], 0):
        I2[i] = tmp_res
        i += 1
    return I2[-1]

# g = 5 / (1 + I2/5)
# for j in range(20):
#     print("total number of segments" + str((1+j) * 10000))
#     print(solve((1+j) * 10000))

segments = 1000
total_length = 10
step_size = total_length/segments
# 1000 segments, define initial conditions
# z2 = np.array([step_size * i for i in range(segments+1)])
I2 = np.zeros([segments+1])
g = np.zeros([segments+1])
I2[0] = 10 ** (-3)
g[0] = 5/(1 + I2[0]/5)
i = 1
for (tmp_res, gain) in ode_solver(dinten_dz, segments, step_size, I2[0], 0):
    I2[i] = tmp_res
    g[i] = gain
    i += 1
print(g.shape)
print(g)
np.savetxt('gain.csv', g, delimiter=',')

# plot the results
# fig, ax1 = plt.subplots()

# plt.xlim(0,10)
# ax1.set_xlabel(r'$z/m$')
# ax1.set_ylabel(r'$Intensity/(W/m^2)$', color='tab:red')
# print(I2[-1:0])
# ax1.plot(z2, I2, 'o',
#             # linestyle='solid',
#             ms=1, color='tab:red',
#             label=r'$\textrm{1000 steps}$')
# # ax1.legend() # loc='center right') 
# ax1.set_yscale("log")

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# color = 'tab:brown'
# # we already handled the x-label with ax1
# ax2.set_ylabel(r"$\textrm{Saturated Gain}/m^{-1}$", color=color)  
# ax2.plot(z2, g, 'o', ms=1, color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# ax2.set_yscale("log")

# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# fig.savefig("q4b.png")