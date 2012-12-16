##############
# test bicycle
##############

import bicycle as bi
from numpy import pi
#from dtk import bicycle as bicycle

#===============
#parameters test
bp = bi.benchmark_parameters() #assertation bicycle.benchmark_parameters
mp = bi.benchmark_to_moore(bp) #assertation bicycle.benchmark_to_moore


#==========
#pitch test
lam = bi.lambda_from_132(mp['rf'], mp['rr'], mp['d1'], mp['d3'], mp['d2'])
#assertation bicycle.lambda_from_abc

lean = 0.0; steer = 0.0
pitch = bi.pitch_from_roll_and_steer(lean, steer, mp['rf'], mp['rr'], 
                                mp['d1'], mp['d2'], mp['d3'])
#assertation bicycle.pitch_from_roll_and_steer(lean, steer, mp['rf'], mp['rr'], 
#                                       mp['d1'], mp['d2'], mp['d3'])


#=========
#basu test

#input
basu_input = bi.basu_table_one_input()
stefen_input = bi.basu_to_stefen_input(basu_input, mp['rr'], bp['lambda'])
#assertation bicycle.basu_to_moore_input

basu_output = bi.basu_table_one_output()
stefen_output = bi.basu_to_stefen_output(basu_output)
