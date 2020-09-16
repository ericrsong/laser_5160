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


def ode_solver(expression, steps, step_size, initial_value, t0, tau=None):
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
        k1 = expression(t0, initial_value, tau)
        k2 = expression(t0 + step_size/2, initial_value + step_size * k1 /2, tau)
        k3 = expression(t0 + step_size/2, initial_value + step_size * k2 /2, tau)
        k4 = expression(t0 + step_size, initial_value + step_size * k3, tau)
        initial_value += step_size * (k1 + 2 * k2 + 2 * k3 + k4) /6
        t0 += step_size
        yield initial_value


def euler(expression, steps, step_size, initial_value, t, tau):
    for i in range(steps):
        k1 = expression(t, initial_value, tau)
        initial_value += step_size * k1
        yield initial_value

def dN1_dt(t, N2, tau2):
    return -N2/tau2 + 1

x = np.array([10** (-3) * i for i in range(10 ** 4)])
y1 = np.zeros([10 ** 4])
y2 = np.zeros([10 ** 4])
i = 0
for tmp_res in ode_solver(dN1_dt, 10** (4), 10** (-3), 0, 0, 1):
    y1[i] = tmp_res
    i += 1

i = 0
for tmp_res in euler(dN1_dt, 10** (4), 10** (-3), 0, 0, 1):
    y2[i] = tmp_res
    i += 1

tau1 = 2
tau2 = 1

t = np.array([0.01 * i for i in range(1200)])

N1 = tau1 * (1 + tau1/(tau2 - tau1) * np.exp(-t/tau1) 
        - tau2/(tau2 - tau1) * np.exp(-t/tau2))
N2 = tau2 * (- np.exp(-t/tau2) + 1)


fig, ax1 = plt.subplots()

plt.xlim(0,30)
ax1.set_xlabel(r'$t/\mu s$')
ax1.set_ylabel(r'$N/R_p N_0$') # , color=color)
ax1.plot(t, N1, linestyle='solid', ms=1, color='tab:green', label=r'$N_1$')
ax1.plot(t, N2, linestyle='solid', ms=1, color='tab:red', label=r'$N_2$')
ax1.plot(x, y1, linestyle='solid', ms=1, color='tab:blue',
            label=r'$\textrm{RK4}$')
ax1.plot(x, y2, linestyle='solid', ms=1, color='b',
            label=r'\textrm{Euler}')


def model(y,t):
    tau = 1
    dydt =  - y/ tau + 1
    return dydt

# initial condition
y0 = 0

# time points
t2 = 5 * np.arange(0,10)

# solve ODE
y = odeint(model,y0,t2)

# plot results
# plt.plot(t,y)
# plt.xlabel('time')
# plt.ylabel('y(t)')
# plt.show()
ax1.plot(t2, y, linestyle='solid', ms=1, color='tab:purple', label=r'$odeint$')

ax1.legend() # loc='center right') 
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig("q2_test_odeint.png")