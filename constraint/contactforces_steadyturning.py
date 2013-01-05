import model as mo
import bicycle as bi
import steadyturning as st

import sympy as sym

from numpy import pi


#=============
#bicycle model
biModel = mo.BicycleModel()

biModel.forcing_full() #forceFull matrix

biModel.contact_forces() #conForceNoncontri

contact_forces = biModel.conForceNoncontri

#states assignment
u1, u3, u6 = biModel.speedsDe
u2, u4, u5 = biModel.speedsInde

T4 = biModel.inputForces[0]
Fx_r, Fy_r, Fx_f, Fy_f = biModel.auxiliaryForces


#===========
# parameters
bp = bi.benchmark_parameters()
mp = bi.benchmark_to_moore(bp)

biModel.parameters_symbols(mp)
para_dict = biModel.parameters


#=============================
# steady turning configuration

#ud: {u1d: 0.0, u2d: 0.0, u3d: 0.0, u4d: 0.0, u5d: 0.0, u6d: 0.0}
ud_dict = st.speeds_zeros(biModel.speedsDerivative)

#u
u_dict = {u2: 0.0, u3: 0.0, u4: 0.0}

class SteadyTurning(object):
    """Steady turning class for equalibrium values and contact forces."""


    def __init__(self, lean, steer):
        """Given lean and steer angles for a steady-turning configuration."""

        #lean = pi/8;  steer = pi/4


        #===================
        # equilibrium values

        #configuration
        print ('configuration')
        q_dict, q_dict_d = st.configuration(lean, steer, mp)
        self.configuration = q_dict
        self.configurationDegree = q_dict_d

        print q_dict, 
        print q_dict_d, '\n'

        #dynamic equations
        print ('dynamic equations')
        dynamic_equ = st.forcing_dynamic_equations(biModel.forceFull, 
                                                    para_dict, q_dict, u_dict)
        self.dynamicEquation = dynamic_equ

        print dynamic_equ, '\n'

        #nonholonomic equations
        print ('nonholonomic equations')
        inde_expression, inde_expression_subs = st.de_by_inde(biModel.nonholonomic, 
                                                        q_dict, para_dict, u_dict)
        self.nonholoEquation = inde_expression_subs

        print inde_expression
        print inde_expression_subs, '\n'

        #combination
        dynamic_nonho_equ = st.dynamic_nonholonomic_equations(inde_expression_subs, 
                                                                dynamic_equ)
        self.dynamicnonholoEquation = dynamic_nonho_equ

        print dynamic_nonho_equ, '\n'

        #Calculations
        print ('Calculation')
        u5_value = sym.solve(dynamic_nonho_equ[0], u5)[0] #choose the negative value

        u_others_dict = {u5: u5_value}

        T4_value = sym.solve(dynamic_nonho_equ[1], T4)[0].subs(u_others_dict)

        u1_value = inde_expression_subs[0].subs(u_others_dict)
        u6_value = inde_expression_subs[2].subs(u_others_dict)

        u_others_dict.update(dict(zip([u1, u6, T4], [u1_value, u6_value, T4_value])))

        self.equilibrium = u_others_dict

        print u_others_dict, '\n'


        #========================================
        # contact forces in each body-fixed coord
        print ('Contact forces')
        contact_forces_st = [value.subs(ud_dict).subs(u_dict).subs(para_dict).subs(q_dict).subs(u_others_dict) 
                          for value in contact_forces]

        self.contactForces = contact_forces_st

        for value in contact_forces_st:
            print value, '\n'

        contact_forces_value = sym.solve(contact_forces_st, [Fx_r, Fy_r, Fx_f, Fy_f])

        self.contactForcesValue = contact_forces_value
