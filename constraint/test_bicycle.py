##############
# test bicycle
##############

import bicycle as bi
#from dtk import bicycle as bicycle

bp = bi.benchmark_parameters() #assertation bicycle.benchmark_parameters
mp = bi.benchmark_to_moore(bp) #assertation bicycle.benchmark_to_moore

lam = bi.lambda_from_132(mp['rf'], mp['rr'], mp['d1'], mp['d3'], mp['d2'])
#assertation bicycle.lambda_from_abc


