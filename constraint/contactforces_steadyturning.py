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



class SteadyTurning(object):
    """Steady turning class for equalibrium values and contact forces."""


    def __init__(self, lean, steer):
        """Given lean and steer angles for a steady-turning configuration."""

        # Bicycle model: biModel
        # forceFull matrix
        # conForceNoncontri
        biModel = mo.BicycleModel()

        self._forceFull = biModel.forcing_full()
        self._contactforcesOrig = biModel.contact_forces()

        # States assignment
        # Parameters
        u1, u3, u6 = biModel._speedsDe
        u2, u4, u5 = biModel._speedsInde
        T4 = biModel._inputForces[0]

        bp = bi.benchmark_parameters()
        mp = bi.benchmark_to_moore(bp)
        self._parameters = mo.strings2symbols(mp, go2type="orsymbols")

        # Steady turning configuration:
        self._ud0s = mo.zeros_dict(biModel._speedsDerivative)
        self._u0s = mo.zeros_dict([u2, u3, u4])

        self._equilibriumSym = [u1, u5, u6, T4]
        self._contactforcesSym = biModel._auxiliaryForces

        # Configuration: e.g. lean = pi/8;  steer = pi/4
        q_dict, q_dict_d = st.configuration(lean, steer, mp)
        self._configuration = q_dict
        self._configurationDegree = q_dict_d

        # Dynamic equations
        # Nonholonomic equations
        dynamic_equ = st.forcing_dynamic_equations(self._forceFull, 
                                            self._parameters, q_dict, self._u0s)
        self._dynamic = dynamic_equ

        inde_expre, inde_expre_subs = st.de_by_inde(biModel._nonholonomic, 
                                            q_dict, self._parameters, self._u0s)
        self._nonho = inde_expre
        self._nonhoSubs = inde_expre_subs

        # Combination (Substitution)
        dynamic_nonho_equ = st.dynamic_nonholonomic_equations(inde_expre_subs, 
                                                                dynamic_equ)
        self._dynamicnonho = dynamic_nonho_equ

    def equi_cal(self):
        """Calculate equilibrium values in the specified steady turning.
        The values encompass u1, u4, u5 and T4."""

        # calculation orders: u5_value, T4_value, u1_value, u6_value
        # Here, should consider more condition, but only choose the negative value
        equi = self._equilibriumSym

        u5_value = sym.solve(self._dynamicnonho[0], equi[1])[0]

        u_others_dict = {equi[1]: u5_value}

        T4_value = sym.solve(self._dynamicnonho[1], equi[3])[0].subs(u_others_dict)

        u1_value = self._nonhoSubs[0].subs(u_others_dict)
        u6_value = self._nonhoSubs[2].subs(u_others_dict)

        u_others_dict.update(dict(zip(equi[:1] + equi[2:], 
                                [u1_value, u6_value, T4_value])))
        self.equilibrium = u_others_dict

        return self.equilibrium

    def contact_force(self):
        """Calculate the contact forces in the defined configuration and 
        equilibrium values. 
        The contact forces are expressed by each body-fixed coordinates."""

        equilibrium = self.equi_cal()

        contact_forces_st = [value.subs(self._ud0s).subs(self._u0s)
                                .subs(self._parameters)
                                .subs(self._configuration).subs(equilibrium) 
                                for value in self._contactforcesOrig]
        self.contactForcesExpre = contact_forces_st

        contact_forces_value = sym.solve(contact_forces_st, self._contactforcesSym)
        self.contactForcesValue = contact_forces_value
