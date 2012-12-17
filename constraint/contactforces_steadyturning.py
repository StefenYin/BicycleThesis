import model as mo
import bicycle as bi
import steadyturning as st
import sympy as sym
from numpy import pi

#=============
#bicycle model
biModel = mo.bicycle_model() #class

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

#ud
ud_dict = st.speeds_zeros(biModel.speedsDerivative)
#ud_dict = {u1d: 0.0, u2d: 0.0, u3d: 0.0, u4d: 0.0, u5d: 0.0, u6d: 0.0}

#u
u_dict = {u2: 0.0, u3: 0.0, u4: 0.0}

#q
lean = pi/8;  steer = pi/4


#===================
# equilibrium values

#configuration
print ('configuration')
q_dict, q_dict_d = st.configuration(lean, steer, mp)
print q_dict, 
print q_dict_d, '\n'


#dynamic equations
print ('dynamic equations')
dynamic_equ = st.forcing_dynamic_equations(biModel.forceFull, para_dict, q_dict, u_dict)
print dynamic_equ, '\n'


#nonholonomic equations
print ('nonholonomic equations')
inde_expression, inde_expression_list = st.de_by_inde(biModel.nonholonomic, 
                                                q_dict, para_dict, u_dict)
print inde_expression
print inde_expression_list, '\n'


#combination
dynamic_nonho_equ = st.dynamic_nonholonomic_equations(inde_expression_list, dynamic_equ)
print dynamic_nonho_equ, '\n'


#Calculations
print ('Calculation')
u5_value = sym.solve(dynamic_nonho_equ[0], u5)[0]  #choose the negative value

u_others_dict = {u5: u5_value}

T4_value = sym.solve(dynamic_nonho_equ[1], T4)[0].subs(u_others_dict)

u1_value = inde_expression_list[0].subs(u_others_dict)
u6_value = inde_expression_list[2].subs(u_others_dict)

u_others_dict.update(dict(zip([u1, u6], [u1_value, u6_value])))

print u_others_dict
print T4_value, '\n'


#========================================
# contact forces in each body-fixed coord
print ('Contact forces')
contact_forces_st = [value.subs(ud_dict).subs(u_dict).subs(para_dict).subs(q_dict).subs(u_others_dict) 
                  for value in contact_forces]

for value in contact_forces_st:
    print value, '\n'

sym.solve(contact_forces_st,[Fx_r, Fy_r, Fx_f, Fy_f])
