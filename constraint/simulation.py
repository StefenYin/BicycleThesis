"""
This module is used for simulating the motion of the bicycle model according
to the ordinary differential equations derived from Kane' method.

The motion will happen in the reference equilibrium values and another steady
state turning equilibrium values in order to obtain the clean states values
from simulation for the calculation of contact forces. Also, the a little
change of initial values will be put in the system for the simulation. Finally,
white noise will be added into the states values for test the robust of the
model and equations of contact forces.

Here, the equations will be calculated from KaneMethod class all from model.py 
module.
"""
from matplotlib.pyplot import (figure, plot, legend, xlabel, ylabel, show,
                                title)

def plotting(q, u, y_ode, note, show=True):
    """
    Parameter
    ---------
    q: list
       a list of states angles.
    u: list
       a list of states rate. 
    y_ode: an array, time by states
       The result from odeint calculation. 
       Rows are time vector, cols are states.
    note: string
       A title description of a plot.
    """

    for i in range(len(q)):
        figure(1)
        plot(t, y_ode[:, i], label='q'+str(i+1))
    for j in range(len(u)):
        figure(2)
        plot(t, y_ode[:, j+len(q)], label='u'+str(j+1))

    figure(1)
    title('Angle performance in %s configuration'%note)
    legend(loc=0)
    xlabel('Time (s)')
    ylabel('Angle (rad)')
        
    figure(2)
    title('Rate performance in %s configuration'%note)
    legend(loc=0)
    xlabel('Time (s)')
    ylabel('Rate (rad/s)')

    if show:
        figure(1).show()
        figure(2).show()

    return figure(1), figure(2)


from scipy.integrate import odeint
from sympy import (Dummy, lambdify, symbols, zeros)
from sympy.physics.mechanics import dynamicsymbols
from numpy.linalg import solve
from numpy import (array, hstack, linspace, pi)

from model import BicycleModel
from bicycle import (benchmark_parameters, benchmark_to_moore, 
                      pitch_from_roll_and_steer)

# Import BicycleModel
biModel = BicycleModel()
left = biModel.mass_matrix_full()
right = biModel.forcing_full()

# symbols
# here, we introduce ua and uad, since forcing equations have them.
# we need to check the forcing?
ua = biModel._auxiliarySpeeds
uad = biModel._auxiliarySpeedsDerivative
ua_zero = {uai:0 for uai in ua}
uad_zero = {uadi:0 for uadi in uad}
q = biModel._coordinatesInde + biModel._coordinatesDe
u = biModel._speedsInde + biModel._speedsDe
T4 = biModel._inputForces[0]
dummy_symbols = [Dummy() for i in q+u]
symbols_dict = dict(zip(q+u, dummy_symbols))
symbols_dict.update({T4: 0.})
v = symbols('v')

# Parameters
bp = benchmark_parameters()
mp = benchmark_to_moore(bp)
para_dict = {}
for key, value in mp.items():
    para_dict.update(dict(zip([symbols(key)], [value])))

# Substitution
MM = left.subs(symbols_dict).subs(para_dict).subs(uad_zero).subs(ua_zero)
For = right.subs(symbols_dict).subs(para_dict).subs(uad_zero).subs(ua_zero)

# Built quick functions for mm and forcing
# Built a function of derivs for q_dot and u_dot
# * here is to convert variables of array (maybe or list) to a tuple for input
# into the lambdify
mm = lambdify(dummy_symbols, MM)
fo = lambdify(dummy_symbols, For)
def derivs(y, t):
    return array(solve(mm(*y), fo(*y))).T[0]

# Initial conditions and time vector defining
# initial condition:
#                  [q1,   q2,   q4,   q3,   u2,   u4,   u5,   u1,   u3,   u6]
#   reference      [0.,   0.,   0.,   lam,  0.,   0.,   v/rr, 0.,   0.,   v/rf]
q1 = 0.; q2 = 0.; q4 = 0.
q3 = pitch_from_roll_and_steer(q2, q4, mp['rf'], mp['rr'], mp['d1'], mp['d2'],
                                mp['d3'])
q_ini = [q1, q2, q4, q3]

v = 5.0 # m/s
u2 = 0.; u4 = 0.; u5 = v/mp['rr']
u1 = 0.; u3 = 0.; u6 = v/mp['rf']
u_ini = [u2, u4, u5, u1, u3, u6]

y0 = array(q_ini + u_ini)
t = linspace(0, 10, 1000)

# Intgration
y = odeint(derivs, y0, t)
y_re = hstack((y[:, :2], y[:,3:4], y[:,2:3], 
               y[:,5:6], y[:,7:9], y[:,4:5], y[:,6:7], y[:,9:]))

# Plotting
figure1, figure2 = plotting(q, u, y_re, 'Reference')
