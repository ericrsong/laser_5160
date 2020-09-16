import numpy as np 

def ode_solver(expression, steps, step_size, initial_value, t0, tau1, tau2):
    """This numerically the ODE using 4th order Runge-Kutta method
        ref: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#

    Parameters:
        expression (func): f(t,y) using the notation in ref above
        steps (int): Total number of steps needed to be integrated
        stepsize (float): The size of each step in t, i.e. Δ
        initial value (float): y(t_0)
        t0 (float): t_0 initial time of the integration
        tau1, tau2 (float): extra parameter for the function expression

    Returns:
        initial value (float): a generator for value of y at each step 
                                of integration
    """
    for i in range(steps):
        k1 = expression(t0, initial_value, tau1, tau2)
        k2 = expression(t0 + step_size/2, 
                        initial_value + step_size * k1 /2, tau1, tau2)
        k3 = expression(t0 + step_size/2,
                        initial_value + step_size * k2 /2, tau1, tau2)
        k4 = expression(t0 + step_size,
                        initial_value + step_size * k3, tau1, tau2)
        initial_value += step_size * (k1 + 2 * k2 + 2 * k3 + k4) /6
        t0 += step_size
        yield initial_value


def dN1_dt(t, N, tau1, tau2):
    """This specifies the ODE that we are solving in this assignment

    Parameters:
        t (float): dummy variable, place-holder for potential future use
        N (nd-array): specifies the occupation number in each state, here it is 
                        denoted in terms of [N1, N2]
        tau1 (float): lifetime of the state |1>
        tau2 (float): lifetime of the state |2>

    Returns:
        (nd-array): approximation of local value of derivative [dN1/dt, dN2/dt]
    """
    return np.array([-N[0]/tau1 + N[1]/tau2, -N[1]/tau2 + 1])

# Input parameters
tau1 = 2
tau2 = 1
init = np.array([0.,0.]) # Initial conditions
total_steps = 10 ** 4
step_size = 10 ** (-3)

# time coordinate of each update
t1 = np.array([step_size * (i+1) for i in range(total_steps)])
y1 = np.zeros([total_steps,2]) # N(t) coordinate for each update [N1(t), N2(t)]
i = 0
for tmp_res in ode_solver(dN1_dt, total_steps, step_size, init, 0, tau1, tau2):
    y1[i] = tmp_res
    i += 1