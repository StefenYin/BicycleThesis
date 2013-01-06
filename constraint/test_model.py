"""
test_model module:
1, importing BicycleModel, benchmark parameters, and zero auxiliary speeds;
2, nonlinear with basu;
3, linearized with benchmark.
"""

import pdb

import bicycle as bi
import model as mo


# Call Class BicycleModel: biModel
# forcing_full, mass_matrix_full
biModel = mo.BicycleModel()

biModel.forcing_full()

biModel.mass_matrix_full()

# Parameters
bp = bi.benchmark_parameters()
mp = bi.benchmark_to_moore(bp)

biModel.parameters_symbols(mp)
para_dict = biModel.parameters

# Zero auxiliary speeds
biModel.auxiliary_speeds_zero()
ua_dict = biModel.auxiliarySpeedsZeros

# Basu-Mandal for nonlinear model
# Input forces or torques: T4 to be zero
# Derivative of zero: deri to be zero since I did not find a good way to wipe it.
# Basu inputs and outputs to stefen ones: input_states_dict, output_dict.
# Calculation: output_cal
# Assertation: assert output_cal == output_dict
T4 = biModel.inputForces[0]
steerTorque = {T4: 0.0}

deri = {'Derivative(0, t)': 0.0}

basu_input = bi.basu_table_one_input()
basu_output = bi.basu_table_one_output()
stefen_input = bi.basu_to_stefen_input(basu_input, mp['rr'], bp['lambda'])
stefen_output = bi.basu_to_stefen_output(basu_output)

biModel.coordinates_dynamicsymbols(stefen_input) 
input_states_dict = biModel.coordinates

biModel.speeds_dynamicsymbols(stefen_output)
output_dict = biModel.speeds

mass_full_nonlin = biModel.mmFull.subs(
                                    ua_dict).subs(
                                    steerTorque).subs(
                                    para_dict).subs(
                                    input_states_dict).subs(
                                    deri)
force_full_nonlin = biModel.forceFull.subs(
                                        ua_dict).subs(
                                        steerTorque).subs(
                                        para_dict).subs(
                                        input_states_dict).subs(
                                        deri)
output_cal = mass_full_nonlin.inv()*force_full_nonlin


# Linearized model
# Configuration
biModel.linearized_reference_configuration(bp['lambda'], mp['rr'], mp['rf'])
linearized_confi = biModel.referenceConfiguration

# Benchmark for Linearization
# mass_full_lin, forcing_lin_A for A matrix
# Assertation: 
# import dtk.bicycle as bicycle
# M, C1, K0, K2 = bicycle.benchmark_par_to_canonical(bp)
# bicycle.benchmark_state_space(M, C1, K0, K2, v, bp['g']) #define v FIRST
mass_full_lin = biModel.mmFull.subs(para_dict).subs(linearized_confi).evalf()

biModel.linearized_a()
forcing_lin_A = biModel.forcingLinA.subs(
                                        ua_dict).subs(
                                        para_dict).subs(
                                        linearized_confi).evalf()
Amat = mass_full_lin.inv() * forcing_lin_A
Am = Amat.extract([1,2,4,5],[1,2,3,4])

# biModel.linearized_b()
# forcing_lin_B = biModel.forcingLinB.subs(
                                        #ua_dict).subs(
                                        #para_dict).subs(
                                        #linearized_confi).evalf()

# Bmat = mass_full_lin.inv() * forcing_lin_B
