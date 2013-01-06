"""
contactforces_steadyturning module:
1, building bicycle model and states, input parameters
2, building a class of SteadyTurning, primarily for contact forces.
"""

import model as mo
import bicycle as bi
import steadyturning as st

import sympy as sym
from numpy import pi


# Bicycle model: biModel
# forceFull matrix
# conForceNoncontri
biModel = mo.BicycleModel()

biModel.forcing_full()

biModel.contact_forces()

contact_forces = biModel.conForceNoncontri

# States assignment
# Parameters
u1, u3, u6 = biModel._speedsDe
u2, u4, u5 = biModel._speedsInde

T4 = biModel._inputForces[0]
Fx_r, Fy_r, Fx_f, Fy_f = biModel._auxiliaryForces

bp = bi.benchmark_parameters()
mp = bi.benchmark_to_moore(bp)

biModel.parameters_symbols(mp)
para_dict = biModel.parameters

# Steady turning configuration:
# ud: {u1d: 0.0, u2d: 0.0, u3d: 0.0, u4d: 0.0, u5d: 0.0, u6d: 0.0}
# u: u2, u3, and u4
ud_dict = st.speeds_zeros(biModel._speedsDerivative)
u_dict = {u2: 0.0, u3: 0.0, u4: 0.0}



class SteadyTurning(object):
    """Steady turning class for equalibrium values and contact forces."""


    def __init__(self, lean, steer):
        """Given lean and steer angles for a steady-turning configuration."""

        # Configuration: e.g. lean = pi/8;  steer = pi/4
        q_dict, q_dict_d = st.configuration(lean, steer, mp)
        self.configuration = q_dict
        self.configurationDegree = q_dict_d

        # Dynamic equations
        # Nonholonomic equations
        dynamic_equ = st.forcing_dynamic_equations(biModel.forceFull, 
                                                    para_dict, q_dict, u_dict)
        self.dynamicEquation = dynamic_equ

        inde_expression, inde_expression_subs = st.de_by_inde(biModel._nonholonomic, 
                                                        q_dict, para_dict, u_dict)
        self.nonholoEquation = inde_expression
        self.nonholoEquationSubs = inde_expression_subs

        # Combination (Substitution)
        dynamic_nonho_equ = st.dynamic_nonholonomic_equations(inde_expression_subs, 
                                                                dynamic_equ)
        self.dynamicnonholoEquation = dynamic_nonho_equ

        # Equilibrium values: u5_value, T4_value, u1_value, u6_value
        # Here, should consider more condition, but only choose the negative value
        u5_value = sym.solve(dynamic_nonho_equ[0], u5)[0]

        u_others_dict = {u5: u5_value}

        T4_value = sym.solve(dynamic_nonho_equ[1], T4)[0].subs(u_others_dict)

        u1_value = inde_expression_subs[0].subs(u_others_dict)
        u6_value = inde_expression_subs[2].subs(u_others_dict)

        u_others_dict.update(dict(zip([u1, u6, T4], [u1_value, u6_value, T4_value])))
        self.equilibrium = u_others_dict

        # Contact forces in each body-fixed coord
        contact_forces_st = [value.subs(ud_dict).subs(u_dict).subs(para_dict).subs(q_dict).subs(u_others_dict) 
                          for value in contact_forces]
        self.contactForces = contact_forces_st

        contact_forces_value = sym.solve(contact_forces_st, [Fx_r, Fy_r, Fx_f, Fy_f])
        self.contactForcesValue = contact_forces_value
