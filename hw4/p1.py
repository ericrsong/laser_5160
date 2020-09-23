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
# ---------------------- default figure setting ends---------------#


def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lc.set_array(np.asarray(c))

    # add lines to axes and rescale
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lc)
    ax.autoscale()
    return lc


def ode_solver_v1(expression_g, expression_I,
                    steps_t, step_size_t, steps_z, step_size_z,
                    tau, g, I_0, I_vec, g_0, I_sat):
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
    for i in range(steps_t): # step refers to the steps in time; not in z
        k1 = expression_g(g, I_vec, tau, I_sat, g_0)
        # t0 += step_size_t
        g += step_size_t * k1
        for index in range(g.shape[0]):
            if g[index] < 0:
                g[index] = 0

        for j in range(steps_z):
            if j == 0:
                k2 = expression_I(g[j], I_0)
                I_vec[j] = I_0 + step_size_z * k2
            else:
                k2 = expression_I(g[j], I_vec[j-1])
                I_vec[j] = I_vec[j-1] +  step_size_z * k2

        yield (g, I_vec)


# differential equation for the intensity
def dI_dz(g, I):
    return g * I


# differential equation for time evolution of saturated gain
def dg_dt(g, I, tau, I_sat, g_0):
    lgth = g.shape
    g_0 = g_0 * np.ones(lgth)
    dgdt = (g_0 - g)/tau - I * g / (I_sat * tau)
    # for index in range(g.shape[0]):
    #     if dgdt[index] < 0:
    #         dgdt[index] = 0
    return  dgdt


steps_t = 80000
step_size_t = 10 ** (-3)
steps_z = 1000
step_size_z = 0.01
tau = 10
g_0 = 5
I_sat = 5
g = g_0 * np.ones([steps_z]) 
I_0 = 10 ** (-3)  #  * np.ones([1000]) 
I_vec = I_0 * np.ones([steps_z]) 

# g2 = np.zeros([steps_t+1, steps_z])
# I2 = np.zeros([steps_t+1, steps_z])
# g2[0] = g
# I2[0] = I_vec
index = 1
for (g,i) in ode_solver_v1(dg_dt, dI_dz,
                                steps_t, step_size_t, steps_z, step_size_z,
                                tau, g, I_0, I_vec, g_0, I_sat):
    output = g
    if index == (80000 - 2):
        print(i)
        print(i[-1])
    # g2[i] = g
    # I2[i] = I_vec
    # print(g)
    index += 1
# print(I2[-1])


######## begin plotting ########################
# # g = 5 / (1 + I2/5)
# xss = np.zeros([steps_z, steps_t+1])
# for i in range(steps_z):
#     xss[i] = np.array([i * 10** (-3) for i in range(10000+1)])

# # xss = np.transpose(xss)
# g2 = np.transpose(g2)
# # print(g2[:,0])
# # print(g2[:,-1])
# index = 0.01 * np.array([i+1 for i in range(1000)])
# fig, ax = plt.subplots()

# xs = 10** (-3) * np.arange(steps_t + 1)
# ax.plot(xs, I2[:,-1])
# ax.set_ylim(0, 1000)
# fig.savefig("output_intensity.png")

# I2 = np.transpose(I2)
# ax.cla()
# begin = 0
# end = 10000
# lc = multiline(xss[begin:end,:], I2[begin:end,:], index[begin:end],
#                 cmap='plasma', lw=2)
# axcb = fig.colorbar(lc)
# axcb.set_label(r"$\textrm{z axis}$")
# ax = plt.gca()
# ax.set_ylim(0, 1000)
# ax.set_xlabel(r'$Time$')
# ax.set_ylabel(r'$Intensity$')
# # ax.set_yscale("log")
# fig.savefig("intensity.png")

# ax.cla()
# begin = 0
# end = 10000
# lc = multiline(xss[begin:end,:], g2[begin:end,:], index[begin:end],
#                 cmap='plasma', lw=2)
# # axcb = fig.colorbar(lc)
# axcb.set_label(r"$\textrm{z axis}$")
# ax = plt.gca()
# ax.set_ylim(0, 6)
# ax.set_xlabel(r'$Time$')
# ax.set_ylabel(r'$Gain$')
# # ax.set_yscale("log")
# fig.savefig("Gain.png")
######## end plotting ########################

# output =g2[:,-1]
print(output.shape)
np.savetxt('gain_time_dependent_time=80ms.csv', output, delimiter=',')