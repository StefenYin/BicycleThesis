##########################
# test steadyturning_funcs
##########################

import bicycle as bi
import steadyturning_funcs as sTurning
import model as mo
import sympy.physics.mechanics as mec

#======================
print ('bicycle model')
#call Class bicycle_model
biModel = mo.bicycle_model()

biModel.forcing_full() #for F_full

#======================================
print ('parameters and states values')

#----------
#parameters
bp = bi.benchmark_parameters()

mp = bi.benchmark_to_moore(bp)

biModel.parameters_symbols(mp)

para_dict = biModel.parameters

#-----------------------------
# test reference configuration
# u2-leanrate, u3-pitchrate, u4-steerrate

u2, u3, u4 = mec.dynamicsymbols('u2 u3 u4')

u_dict = sTurning.speeds_zeros([u2, u4])

lean = 0.0; steer = 0.0


#========================
print ('steady turning')

#configuration
q_dict, q_dict_d = sTurning.configuration(lean, steer, mp)
print q_dict, 
print q_dict_d, '\n'

#dynamic equations
dynamic_equ = sTurning.forcing_dynamic_equations(biModel.forceFull, para_dict, q_dict, u_dict)
print dynamic_equ, '\n'

#nonholonomic equations
inde_expression, inde_expression_list = sTurning.de_by_inde(biModel.nonholonomic, 
                                                q_dict, para_dict, u_dict)
print inde_expression
print inde_expression_list, '\n'

#substitution
dynamic_nonho_equ = sTurning.dynamic_nonholonomic_equations(inde_expression_list, dynamic_equ)
print dynamic_nonho_equ, '\n'
