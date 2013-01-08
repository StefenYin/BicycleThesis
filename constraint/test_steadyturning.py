"""
test_steadyturning module:
1, building BicycleModel and importing parameters;
2, trying steady turning functions in reference configuration.
"""

from mo import (BicycleModel, strings2symbols)
from bicycle import (benchmark_parameters, benchmark_to_moore)
from steadyturning import (configuration, forcing_dynamic_equations,
                            de_by_inde, dynamic_nonholonomic_equations
                            )


# Call Class BicycleModel: biModel
# forcing_full
biModel = BicycleModel()

forceFull = biModel.forcing_full()

# Parameters and states values:
# Parameters
bp = benchmark_parameters()
mp = benchmark_to_moore(bp)
para_dict = strings2symbols(mp, go2type="orsymbols")

# Test reference configuration
# u2-leanrate, u3-pitchrate, u4-steerrate
u_dict = strings2symbols(biModel._speeds[1:4], go2type="dysymbols")

lean = 0.0; steer = 0.0

# Steady turning:
# Configuration
# Dynamic equations: dynamic_equ
# Nonholonomic equations: inde_expression, inde_expression_subs
# Substitution: dynamic_nonho_equ
q_dict, q_dict_d = configuration(lean, steer, mp)

dynamic_equ = forcing_dynamic_equations(forceFull, 
                                                para_dict, q_dict, u_dict)

inde_expression, inde_expression_subs = de_by_inde(biModel._nonholonomic, 
                                                q_dict, para_dict, u_dict)

dynamic_nonho_equ = dynamic_nonholonomic_equations(inde_expression_subs, dynamic_equ)
