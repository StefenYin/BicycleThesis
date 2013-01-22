"""
This module is used for simulating the motion of the bicycle model according
to the ordinary differential equations derived from Kane' method.

The motion will happen in the reference equilibrium values and another steady
state turning equilibrium values in order to obtain the clean states values
from simulation for the calculation of contact forces. Also, the a little
change of initial values will be put in the system for the simulation. Finally,
white noise will be added into the states values for test the robust of the
model and equations of contact forces.

Here, the equations will be calculated from KaneMethod class all from 
contactforces_steadyturning.py module.
"""
