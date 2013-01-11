"""
test_model module:
1, importing BicycleModel, benchmark parameters, and zero auxiliary speeds;
2, nonlinear with basu;
3, linearized with benchmark.
"""

from model import (BicycleModel, strings2symbols, zeros_dict)
from bicycle import (benchmark_parameters, benchmark_to_moore,
                    basu_table_one_input, basu_table_one_output,
                    basu_to_stefen_input, basu_to_stefen_output
                    )

# Call Class BicycleModel: biModel
# forcing_full, mass_matrix_full, para_dict, ua_dict
# Go to parameters
biModel = BicycleModel()

forceFull = biModel.forcing_full()
mmFull = biModel.mass_matrix_full()

bp = benchmark_parameters()
mp = benchmark_to_moore(bp)
para_dict = strings2symbols(mp, go2type="orsymbols")
ua_dict = zeros_dict(biModel._auxiliarySpeeds)

# Basu-Mandal for nonlinear model
# Input forces or torques: T4 to be zero
# Derivative of zero: deri to be zero since I did not find a good way to wipe it.
# Basu inputs and outputs to stefen ones: input_states_dict, output_dict.
# Calculation: output_cal
# Assertation: assert output_cal == output_dict
T4 = biModel._inputForces[0]
steerTorque = {T4: 0.}

deri = {'Derivative(0, t)': 0.}

basu_input = basu_table_one_input()
basu_output = basu_table_one_output()
stefen_input = basu_to_stefen_input(basu_input, mp['rr'], bp['lambda'])
stefen_output = basu_to_stefen_output(basu_output)

input_states_dict = strings2symbols(stefen_input, go2type="dysymbols") 
output_dict = strings2symbols(stefen_output, go2type="dysymbols")

mass_full_nonlin = mmFull.subs(
                                ua_dict).subs(
                                steerTorque).subs(
                                para_dict).subs(
                                input_states_dict).subs(
                                deri)
force_full_nonlin = forceFull.subs(
                                    ua_dict).subs(
                                    steerTorque).subs(
                                    para_dict).subs(
                                    input_states_dict).subs(
                                    deri)
output_cal = mass_full_nonlin.inv()*force_full_nonlin


# Linearized model
# Configuration
linearized_confi = biModel.linearized_reference_configuration(bp['lambda'], 
                                                        mp['rr'], mp['rf'])

# Benchmark for Linearization
# mass_full_lin, forcing_lin_A for A matrix
# Assertation: 
# import dtk.bicycle as bicycle
# M, C1, K0, K2 = bicycle.benchmark_par_to_canonical(bp)
# bicycle.benchmark_state_space(M, C1, K0, K2, v, bp['g']) #define v FIRST
mass_full_lin = mmFull.subs(para_dict).subs(linearized_confi).evalf()

forcingLinA = biModel.linearized_a()
forcing_lin_A = forcingLinA.subs(
                                ua_dict).subs(
                                para_dict).subs(
                                linearized_confi).evalf()
Amat = mass_full_lin.inv() * forcing_lin_A
Am = Amat.extract([1,2,4,5],[1,2,3,4])
