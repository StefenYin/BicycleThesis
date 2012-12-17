#######################################
# test model, nonlinear with basu
#             linearized with benchmark
########################################

import bicycle as bi
import model as mo
import sympy.physics.mechanics as mec
import pdb

#======================
print ('bicycle model')
#call Class bicycle_model
biModel = mo.bicycle_model()

biModel.forcing_full() #for F_full

biModel.mass_matrix_full() #for MM_full


#====================
print ('parameters')

bp = bi.benchmark_parameters()

mp = bi.benchmark_to_moore(bp)

biModel.parameters_symbols(mp)

para_dict = biModel.parameters


#===============================
print ('zero auxiliary speeds')
biModel.auxiliary_speeds_zero()

ua_dict = biModel.auxiliarySpeedsZeros



#==================================
print ('basu for nonlinear model')

#steer torque
T4 = mec.dynamicsymbols('T4')
steerTorque = {T4: 0.0}

#derivative of zero
deri = {'Derivative(0, t)': 0.0}

#basu inputs and outputs to stefen ones
basu_input = bi.basu_table_one_input()
basu_output = bi.basu_table_one_output()

stefen_input = bi.basu_to_stefen_input(basu_input, mp['rr'], bp['lambda'])
stefen_output = bi.basu_to_stefen_output(basu_output)

biModel.coordinates_dynamicsymbols(stefen_input) #actually, stefen_input includes 
                                            #coordinates and speeds of basu inputs
input_states_dict = biModel.coordinates

#for check
biModel.speeds_dynamicsymbols(stefen_output) #actually, stefen_output includes
                                            #speeds differentiation
output_dict = biModel.speeds

#for calculation of output
mass_full_nonlin = biModel.mmFull.subs(ua_dict).subs(steerTorque).subs(para_dict).subs(input_states_dict).subs(deri)
force_full_nonlin = biModel.forceFull.subs(ua_dict).subs(steerTorque).subs(para_dict).subs(input_states_dict).subs(deri)

output_cal = mass_full_nonlin.inv()*force_full_nonlin



#==========================
print ('linearized model')

#configuration
biModel.linearized_reference_configuration(bp['lambda'], mp['rr'], mp['rf'])

linearized_confi = biModel.referenceConfiguration

#linearization
mass_full_lin = biModel.mmFull.subs(para_dict).subs(linearized_confi).evalf()

biModel.linearized_a()
#biModel.linearized_b()

forcing_lin_A = biModel.forcingLinA.subs(ua_dict).subs(para_dict).subs(linearized_confi).evalf()
#forcing_lin_B = biModel.forcingLinB.subs(ua_dict).subs(para_dict).subs(linearized_confi).evalf()

Amat = mass_full_lin.inv() * forcing_lin_A
#Bmat = mass_full_lin.inv() * forcing_lin_B

Am = Amat.extract([1,2,4,5],[1,2,3,4])
