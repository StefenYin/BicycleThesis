"""
contactforces_steadyturning module:
1, Configuration & Parameters -> turning radius, total center of mass, relative
   position;
2, & dynamic equations, nonholonomic equations -> equilibrium values, 
   contact forces.
3, Equilibrium values checks, encompassing 1) turning radius check; 2) rear
   wheel rate unsolvable; 3) rear wheel rate solved to have imaginary values.
   4) stability of equilibrium values in the defined configuration.
"""

import model as mo
import bicycle as bi
import steadyturning as st
import geometry as ge

import sympy as sym
from numpy import (pi, array, zeros, sqrt)



class SteadyTurning(object):
    """Steady turning class for equalibrium values and contact forces."""


    def __init__(self, lean, steer, bicycleParameters = None):
        """Given lean and steer angles for a steady-turning configuration.

        Parameter
        ---------
        lean, steer: float
            Angles for defining a steady turning configuration.
        bicycleParameters: a dictionary
            The parameters is the dictionary is expressed in Moore set, not 
            benchmark set.
            Default parameters are benchmark parameters.

        """

        # Bicycle model: biModel
        # forceFull matrix
        # conForceNoncontri
        # contact points relative position in a list, but from ways to obtain it
        # individual centers of mass of four rigid bodies
        biModel = mo.BicycleModel()
        self._biModel = biModel

        self._mmFull = biModel.mass_matrix_full()
        self._forcingLinA = biModel.linearized_a()
        self._forceFull = biModel.forcing_full()
        self._contactforcesOrig = biModel.contact_forces()
        self._contactPosition = (biModel._contact_posi_pq + 
                                    biModel._contact_posi_dfn)
        self._turningRadiusSym = biModel._turningRadiusSym
        self._bodies_dn_A = biModel._bodies_dn_A
        self._massSym = biModel._massSym
        self._contactforcesSym = biModel._auxiliaryForces

        # States assignment
        # Parameters
        u1, u3, u6 = biModel._speedsDe
        u2, u4, u5 = biModel._speedsInde
        T4 = biModel._inputForces[0]

        if bicycleParameters is None:
            bp = bi.benchmark_parameters()
            mp = bi.benchmark_to_moore(bp)
        else:
            mp = bicycleParameters
        self._parameters = mo.strings2symbols(mp, go2type="orsymbols")

        # Steady turning configuration:
        # Configuration: e.g. lean = pi/8;  steer = pi/4
        self._ud0s = mo.zeros_dict(biModel._speedsDerivative)
        self._u0s = mo.zeros_dict([u2, u3, u4])
        self._equilibriumSym = [u1, u5, u6, T4]

        q_dict, q_dict_d = st.configuration(lean, steer, mp)
        self._configuration = q_dict
        self._configurationDegree = q_dict_d

        # Total center of mass(including relative position of fn to dn)
        # Turning radius
        self.total_com()
        self.turning_radius()

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

        # Equilibrium values
        self.equi_cal()

    def equi_cal(self):
        """Calculate equilibrium values in the specified steady turning.
        The values encompass u1, u4, u5 and T4."""

        # calculation orders: u5_value, T4_value, u1_value, u6_value
        # Here, should consider more condition, but only choose the negative value
        u1, u5, u6, T4 = self._equilibriumSym

        u5sqr_value = sym.solve(self._dynamicnonho[0], u5**2)
        if u5sqr_value == []:
            raise ValueError("\nOops! Rear wheel rate in configuration {0} \
cannot be solved. Please select valid configuration according to the plot \
from <General steady turning of a benchmark bicycle model> by Luke. \
Good luck!\n".format(self._configuration))
        elif u5sqr_value[0] < 0:
            raise ValueError("Oops! The steady turning in your configuration \
{0} seems Infeasible since no real value appears in rear wheel rate. Please \
check the B8_num_st and try another valid configuraiton according to the plot \
from <General steady turning of a benchmark bicycle model> by Luke.\n"\
.format(self._configuration))
        else:
            u5_value = -sym.sqrt(u5sqr_value[0])
            print ("\nIt passed feasible steady turning checking and already \
solved the equilibrium values, but it still needs to be checked by eigenvalues \
of the configuration...\n")

        u_others_dict = {u5: u5_value}

        T4_value = sym.solve(self._dynamicnonho[1], T4)[0].subs(u_others_dict)

        u1_value = self._nonhoSubs[0].subs(u_others_dict)
        u6_value = self._nonhoSubs[2].subs(u_others_dict)

        u_others_dict.update(dict(zip([u1, u6, T4], 
                                      [u1_value, u6_value, T4_value])))
        self._equilibrium_u = u_others_dict
        self._equilibrium_u.update(self._u0s)


    def contact_force(self):
        """Calculate the contact forces in the defined configuration and 
        equilibrium values. 
        The contact forces are expressed by each body-fixed coordinates."""

        equilibrium = self._equilibrium_u

        contact_forces_st = [value.subs(self._ud0s).subs(self._u0s)
                                .subs(self._parameters)
                                .subs(self._configuration).subs(equilibrium).doit() 
                                for value in self._contactforcesOrig]
        self.contactForcesExpre = contact_forces_st

        contact_forces_value = sym.solve(contact_forces_st, self._contactforcesSym)
        self.contactForcesValue = contact_forces_value

        return self.contactForcesValue


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
        self._totalComA123 = cT.tolist()


    def turning_radius(self):
        """Returns the turning radius of two wheels and total center of mass
        after defining a steady turning configuration."""

        contact_position = [value.subs(self._parameters)
                            .subs(self._configuration) 
                            for value in self._contactPosition]
        self._fn_dn_A1 = contact_position[2]
        self._fn_dn_A2 = contact_position[3]

        Rr, Rf = self._turningRadiusSym

        turn_radius = sym.solve([contact_position[0] - contact_position[2], 
                                contact_position[1] - contact_position[3]],
                                self._turningRadiusSym)
        
        if turn_radius == []:
            print ("It seems the configuration {0} that you are building is \
not going to generate a steady turning. Maybe you need to try another valid \
configuration.\n".format(self._configuration))
            pass
        else:
            print ("It seems the configuration {0} that you are building is \
right. At lease it has turning radius at this point, but it is still being \
checked...\n".format(self._configuration))

            self._turningRadiusRearGeo = turn_radius[Rr]
            self._turningRadiusFrontGeo = turn_radius[Rf]
            self._turningRadiusCom = sqrt(float(self._totalComA123[0]**2 + 
                                    (turn_radius[Rr] - self._totalComA123[1])**2))


    def eigenvalues(self):
        """Returns the eigenvalues of state A matrix in the defined steady
        turning configuration."""

        self._equilibrium_qu = {}
        self._equilibrium_qu.update(self._configuration)
        self._equilibrium_qu.update(self._equilibrium_u)

        mass = self._mmFull.subs(self._parameters).subs(self._equilibrium_qu)
        lin_forcing_a = self._forcingLinA.subs(self._parameters).subs(self._equilibrium_qu)
        self._AmatFull = mass.inv() * lin_forcing_a
        self._Amat = self._AmatFull.extract([4,5,6], [3,4,5])
        eig1, eig2, eig3 = self._Amat.eigenvals(rational = False)

        self._eigs = [complex(eigvalue) for eigvalue in [eig1, eig2, eig3]]

        ei_signs = [eigs.real <= 0. for eigs in self._eigs]

        if False not in ei_signs:
            print ("It passed the eigenvalues checking. You can see the \
eigenvalues in ._eigs. Overall, the configuration {0} seems already generates \
a steady turning.\n".format(self._configuration))
        else:
            print ("Some eigenvalues are positive which means the steady \
turning in the configuration {0} is not stable or kinematically infeasible. \
Please check the eigenvalue and try another configuration.\n".format(self._configuration))
