"""
test_steadyturning module:
1, building BicycleModel and importing parameters;
2, trying steady turning functions in reference configuration.
"""

import bicycle as bi
import steadyturning as sTurning
import model as mo

import sympy.physics.mechanics as mec


# Call Class BicycleModel: biModel
# forcing_full
biModel = mo.BicycleModel()

forceFull = biModel.forcing_full()

# Parameters and states values:
# Parameters
bp = bi.benchmark_parameters()
mp = bi.benchmark_to_moore(bp)
para_dict = mo.strings2symbols(mp, go2type="orsymbols")

# Test reference configuration
# u2-leanrate, u3-pitchrate, u4-steerrate
u2, u3, u4 = mec.dynamicsymbols('u2 u3 u4')

u_dict = sTurning.speeds_zeros([u2, u4])

lean = 0.0; steer = 0.0

# Steady turning:
# Configuration
# Dynamic equations: dynamic_equ
# Nonholonomic equations: inde_expression, inde_expression_subs
# Substitution: dynamic_nonho_equ
q_dict, q_dict_d = sTurning.configuration(lean, steer, mp)

dynamic_equ = sTurning.forcing_dynamic_equations(forceFull, 
                                                para_dict, q_dict, u_dict)

inde_expression, inde_expression_subs = sTurning.de_by_inde(biModel._nonholonomic, 
                                                q_dict, para_dict, u_dict)

dynamic_nonho_equ = sTurning.dynamic_nonholonomic_equations(inde_expression_subs, dynamic_equ)
