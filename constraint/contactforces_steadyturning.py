"""
contactforces_steadyturning module:
1, building bicycle model and states, input parameters
2, building a class of SteadyTurning, primarily for contact forces.
"""

import model as mo
import bicycle as bi
import steadyturning as st
import geometry as ge

import sympy as sym
from numpy import (pi, array, zeros)



class SteadyTurning(object):
    """Steady turning class for equalibrium values and contact forces."""


    def __init__(self, lean, steer):
        """Given lean and steer angles for a steady-turning configuration."""

        # Bicycle model: biModel
        # forceFull matrix
        # conForceNoncontri
        # contact points relative position in a list, but from ways to obtain it
        # individual centers of mass of four rigid bodies
        biModel = mo.BicycleModel()

        self._forceFull = biModel.forcing_full()
        self._contactforcesOrig = biModel.contact_forces()
        self._contactPosition = (biModel._contact_posi_pq + 
                                    biModel._contact_posi_dfn)
        self._turningRadiusSym = biModel._turningRadiusSym
        self._bodies_dn_A = biModel._bodies_dn_A
        self._massSym = biModel._massSym

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

        # Turning radius
        self.turning_radius()
        self.total_com()

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

        return self.contactForcesValue


    def turning_radius(self):
        """Returns the turning radius of two wheels after defining a steady
        turning configuration."""

        contact_position = [value.subs(self._parameters)
                            .subs(self._configuration) 
                            for value in self._contactPosition]
        self._fn_dn_A1 = contact_position[2]
        self._fn_dn_A2 = contact_position[3]

        Rr, Rf = self._turningRadiusSym

        turn_radius = sym.solve([contact_position[0] - contact_position[2], 
                                contact_position[1] - contact_position[3]],
                                self._turningRadiusSym)

        self._turningRadiusRearGeo = turn_radius[Rr]
        self._turningRadiusFrontGeo = turn_radius[Rf]


    def total_com(self):
        """Returns total center of mass, position and mass, of four rigid
        bodies of a bicycle."""

        mp = self._parameters
        mc, md, me, mf = self._massSym

        coordinates = self._bodies_dn_A
        row, col = coordinates.shape
        coordinates_value = zeros((row, col))

        for i in range(row):
            for j in range(col):
                value = coordinates[i][j].subs(mp).subs(self._configuration)
                coordinates_value[i][j] = value

        masses = array([[mp[mc], mp[md], mp[me], mp[mf]]])
        mT, cT = ge.total_com(coordinates_value, masses)

        self._totalMass = mT
        self._totalcomA123 = cT.tolist()
