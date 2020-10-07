import numpy as np 
import sympy as sp
from sympy.abc import R, L 
import matplotlib.pyplot as plt
import matplotlib as plot
plot.use('Agg')
# from scipy.integrate import odeint

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

def free_space(L):
    return sp.Matrix([[1, L],[0, 1]])


def parallel_plate(L, n):
    return sp.Matrix([[1, L/n],[0, 1]])


def mirror(R):
    return sp.Matrix([[1, 0],[-2/R, 1]])

# define the total ABCD matrix without the rod
def round_trip(free_space, mirror, R, L):
    return free_space(L)*mirror(R)*free_space(L)

# define the total ABCD matrix with the rod
def round_trip_with_rod(free_space, mirror, R, L):
    ABCD =  (free_space(0.1) * parallel_plate(1,1.82) * free_space(L - 1.1) * 
            mirror(R) * free_space(L - 1.1) * parallel_plate(1,1.82) * 
            free_space(0.1))
    return ABCD

# check if the stability condition is satisfied
def if_stable(RT, L, value=10):
    A_D = RT[0, 0].subs(L, value) + RT[1, 1].subs(L, value)
    if A_D < 2 and A_D > -2:
        return True
    elif A_D == 2 or A_D == -2:
        print("The system is critically stable")
        return True
    else:
        return False


def r_waist(RT, L, value):
    if if_stable(RT, L,value):
        # wavelength of light in cavity
        lambda_1 = 1.064 * 10 ** (-4)
        # trace of the ABCD matrix
        A_D = RT[0, 0].subs(L, value) + RT[1, 1].subs(L, value)
        abs_C = sp.Abs(RT[1,0].subs(L,value))
        # find the expression for omega_0
        expres = sp.sqrt(lambda_1 * sp.sqrt(4- A_D ** 2)/(2 * np.pi * abs_C))
        return expres
    else:
        raise ValueError("The system is not stable for input of value L= " +
                            str(value))


R = 10
RT = round_trip_with_rod(free_space, mirror, R, L)

x = np.array([0.001 * i for i in range(1101,12000)])
xcoor = []
ycoor = []
for i in range(x.shape[0]):
    try:
        ycoor.append(r_waist(RT, L, x[i]))
        xcoor.append(x[i])
        # print(ycoor)
    except:
        print("The system is not stable for input of value L= " +
                            str(x[i]))
# plot the results
fig, ax = plt.subplots()
ax.plot(xcoor, np.array(ycoor)*10)
ax.set_xlabel(r'$L/cm$')
ax.set_ylabel(r'$Spot Size/mm$')
fig.savefig("spot_size.png")

# def r_waist_sym(RT, L, value):
#     if if_stable(RT, L,value):
#         lambda_1 = 1.064 * 10 ** (-4)
#         A_D = RT[0, 0] + RT[1, 1]
#         abs_C = sp.Abs(RT[1,0])
#         expres = sp.sqrt(lambda_1 * sp.sqrt(4- A_D ** 2)/(2 * np.pi * abs_C))
#         return expres
#     else:
#         raise ValueError("The system is not stable for input of value L= " +
#                             str(value))
# print("the waist is " + str(r_waist(RT, L, 5.)) + "cm")
# print(r_waist_sym(RT, L, 10))
