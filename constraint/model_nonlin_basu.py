################################################
# nonlinear constraint benchmark bicycle model
# check with basu-madal nonlinear values
#################################################

import sympy as sym
from numpy import *
import sympy.physics.mechanics as mec

import os

from dtk import bicycle as bi
mec.Vector.simp = False
mec.mechanics_printing()

##################
# Reference Frames
##################

# Newtonian Frame
N = mec.ReferenceFrame('N', indices=('1', '2', '3'))
       
# Yaw Frame
A = mec.ReferenceFrame('A', indices=('1', '2', '3'))
        
# Roll Frame
B = mec.ReferenceFrame('B', indices=('1', '2', '3'))
        
# Pitch & Bicycle Frame
C = mec.ReferenceFrame('C', indices=('1', '2', '3'))

# Steer & Fork/Handlebar Frame
E = mec.ReferenceFrame('E', indices=('1', '2', '3'))
       
       
# Rear Wheel Frame
D = mec.ReferenceFrame('D', indices=('1', '2', '3'))
       
# Front Wheel Frame
F = mec.ReferenceFrame('F', indices=('1', '2', '3'))
        

####################################
# Generalized Coordinates and Speeds
####################################

# q1,u1: frame yaw
# q2,u2: frame roll
# q3,u3: frame pitch 
# q4,u4: steering rotation
# u5: rear wheel ang. vel.
# u6: front wheel ang. vel.
q1, q2, q3, q4 = mec.dynamicsymbols('q1 q2 q3 q4')
q1d, q2d, q3d, q4d = mec.dynamicsymbols('q1 q2 q3 q4', 1)

u1, u2, u3, u4 = mec.dynamicsymbols('u1 u2 u3 u4')
u5, u6 = mec.dynamicsymbols('u5 u6')

u1d, u2d, u3d, u4d = mec.dynamicsymbols('u1 u2 u3 u4', 1)
u5d, u6d = mec.dynamicsymbols('u5 u6', 1)


#################################
# Orientation of Reference Frames
#################################


# bicycle frame yaw
A.orient(N, 'Axis', [q1, N['3']])
# bicycle frame roll
B.orient(A, 'Axis', [q2, A['1']])
# bicycle frame pitch
C.orient(B, 'Axis', [q3, B['2']])
# fork/handlebar steer
E.orient(C, 'Axis', [q4, C['3']])

###########
# Constants
###########

# geometry
# rF: radius of front wheel
# rR: radius of rear wheel
# d1: the perpendicular distance from the steer axis to the center
#     of the rear wheel (rear offset)
# d2: the distance between wheels along the steer axis
# d3: the perpendicular distance from the steer axis to the center
#     of the front wheel (fork offset)
# l1: the distance in the c1> direction from the center of the rear
#     wheel to the frame center of mass
# l2: the distance in the c3> direction from the center of the rear
#     wheel to the frame center of mass
# l3: the distance in the e1> direction from the front wheel center to
#     the center of mass of the fork
# l4: the distance in the e3> direction from the front wheel center to
#     the center of mass of the fork

rF, rR = sym.symbols('rF rR')
d1, d2, d3 = sym.symbols('d1 d2 d3')
l1, l2, l3, l4 = sym.symbols('l1 l2 l3 l4')

# acceleration due to gravity
g = sym.symbols('g')

# mass
mc, md, me, mf = sym.symbols('mc md me mf')

# inertia
ic11, ic22, ic33, ic31 = sym.symbols('ic11 ic22 ic33 ic31') #rear frame
id11, id22 = sym.symbols('id11 id22')  #rear wheel
ie11, ie22, ie33, ie31 = sym.symbols('ie11 ie22 ie33 ie31')  #front frame
if11, if22 = sym.symbols('if11 if22') #front wheel

###########
# Specified
###########

# control torques
#T2--roll tor, T4--steer tor, T5--rear wheel tor
T4 = mec.dynamicsymbols('T4')

####################################       
# Kinematical Differential Equations
####################################

kinematical = [q1d - u1,
               q2d - u2,
               q3d - u3,
               q4d - u4]
 
####################
# Angular Velocities
####################

A.set_ang_vel(N, u1 * N['3'])
B.set_ang_vel(A, u2 * A['1'])
C.set_ang_vel(B, u3 * B['2'])
E.set_ang_vel(C, u4 * C['3'])

D.set_ang_vel(C, u5 * C['2'])
F.set_ang_vel(E, u6 * E['2'])
 

print('ready for special unit vectors; and points and velocities')

######################
#special unit vectors
######################
#pitch back for the unit vector g_3 along front wheel radius;

g_3 =  (mec.express(A['3'], E) - mec.dot(E['2'], A['3'])*E['2']).normalize() 
#another way: g_3 = E['2'].cross(A['3']).cross(E['2']).normalize() 


#roll back for longitudinal and lateral unit vector of front wheel

long_v = mec.cross (E['2'], A['3']).normalize()
lateral_v = mec.cross (A['3'], long_v).normalize() 


#########################
#points and velocities
#########################

####################rear wheel contact point, dn
dn = mec.Point('dn')
dn.set_vel(N,0.0)

# rear wheel center, do
do = dn.locatenew('do', -rR * B['3']) 
#do = dn.locatenew('do', -rtr * A['3'] - rR * B['3']) # rear wheel center with rear tire radius rtr

do.v2pt_theory(dn, N, D)

t = sym.symbols('t')
do.set_acc(N, do.vel(N).diff(t, B) + mec.cross(B.ang_vel_in(N), do.vel(N))) 


# rear frame center
co = mec.Point('co')
co.set_pos(do, l1 * C['1'] + l2 * C['3'])

co.v2pt_theory(do, N, C)
co.a2pt_theory(do, N, C)

# steer axis point
ce = mec.Point('ce')
ce.set_pos(do, d1 * C['1'])

ce.v2pt_theory(do, N, C)
ce.a2pt_theory(do, N, C)


####################Front wheel contact point, fn
fn = mec.Point('fn')
fn.set_vel(N, 0.0)

#front wheel center
fo = fn.locatenew('fo', -rF * g_3)
#fo = fn.locatenew('fo', -ftr * A['3'] - rF * g_3) # rear wheel center with rear tire radius rtr.

fo.v2pt_theory(fn, N, F)

fo.set_acc(N, fo.vel(N).diff(t, E) + mec.cross(E.ang_vel_in(N), fo.vel(N)))


# front frame center
eo = mec.Point('eo')
eo.set_pos(fo, l3 * E['1'] + l4 * E['3'])

eo.v2pt_theory(fo, N, E)
eo.a2pt_theory(fo, N, E)

#steer axis foot; one from rear wheel center, the other from front wheel center
SAF = do.locatenew('SAF', d1 * C['1'] + d2 * E['3'])
SAF.set_pos(fo, -d3 * E['1'])

#velociy of SAF in two ways
#v_SAF_1 = v2pt_theory(do, N, C)
#v_SAF_2 = v2pt_theory(fo, N, E)
v_SAF_1 = do.vel(N) + mec.cross(C.ang_vel_in(N), SAF.pos_from(do))
v_SAF_2 = fo.vel(N) + mec.cross(E.ang_vel_in(N), SAF.pos_from(fo))

print('ready for constraints; inertia; rigid bodies; bodylist')
####################
# Motion Constraints
####################
holonomic = [fn.pos_from(dn).dot(A['3'])]

nonholonomic = [(v_SAF_1-v_SAF_2).dot(uv) for uv in E]
                                        
#########
# Inertia
#########

Ic = mec.inertia(C, ic11, ic22, ic33, 0.0, 0.0, ic31)
Id = mec.inertia(C, id11, id22, id11, 0.0, 0.0, 0.0) #rear wheel
Ie = mec.inertia(E, ie11, ie22, ie33, 0.0, 0.0, ie31)
If = mec.inertia(E, if11, if22, if11, 0.0, 0.0, 0.0) #front wheel


##############
# Rigid Bodies
##############
rearFrame_inertia = (Ic, co)
rearFrame=mec.RigidBody('rearFrame',co,C,mc,rearFrame_inertia)

rearWheel_inertia = (Id, do)
rearWheel=mec.RigidBody('rearWheel',do,D,md,rearWheel_inertia)

frontFrame_inertia = (Ie, eo)
frontFrame=mec.RigidBody('frontFrame',eo,E,me,frontFrame_inertia)

frontWheel_inertia = (If, fo)
frontWheel=mec.RigidBody('frontWheel',fo,F,mf,frontWheel_inertia)

bodyList = [rearFrame, rearWheel, frontFrame, frontWheel]


###########################
# Generalized Active Forces
###########################

# gravity
Fco = (co, mc * g * A['3'])
Fdo = (do, md * g * A['3'])
Feo = (eo, me * g * A['3'])
Ffo = (fo, mf * g * A['3'])


forceList = [Fco, Fdo, Feo, Ffo]

###############
# Kane's Method
###############

kane = mec.KanesMethod(N, q_ind=[q1, q2, q4], u_ind=[u2, u4, u5], kd_eqs=kinematical, q_dependent=[q3], configuration_constraints=holonomic, u_dependent=[u1, u3, u6], velocity_constraints=nonholonomic)
#reminder: u1--yaw rate, u2--roll rate, u3--pitch rate, u4--steer rate, u5--rear wheel ang. vel., u6--front wheel ang. vel.

#kane.kindiffeq(kinematical)
(fr, frstar)= kane.kanes_equations(forceList, bodyList)
kdd = kane.kindiffdict()

bp = bi.benchmark_parameters()

mp = bi.benchmark_to_moore(bp)

rr = bp['rR']
lam = bp['lam']

basu_input = bi.basu_table_one_input()

mo_in = bi.basu_to_moore_input(basu_input, rr, lam)

#not include the l1, l2, mc and ic because of the rider 
#which will be added into instrumented bicycle later.
para_dict={l1: mp['l1'], l2: mp['l2'], l3: mp['l3'], l4:mp['l4'], 
        mc: mp['mc'], md:mp['md'], me:mp['me'], mf:mp['mf'], rR:mp['rr'], rF:mp['rf'], 
        g:mp['g'], d1: mp['d1'], d2: mp['d2'], d3: mp['d3'],
        id11: mp['id11'], id22: mp['id22'], if11: mp['if11'], if22: mp['if22'],
        ic11: mp['ic11'], ic22: mp['ic22'], ic33: mp['ic33'], ic31: mp['ic31'],
        ie11: mp['ie11'], ie22: mp['ie22'], ie33: mp['ie33'], ie31: mp['ie31']}

in_dict={q1: mo_in['q3'], q2: mo_in['q4'], q3: mo_in['q5'], q4: mo_in['q7'],
        u1: mo_in['u3'], u2: mo_in['u4'], u3: mo_in['u5'], u4: mo_in['u7'],
        u5: mo_in['u6'], u6: mo_in['u8']}


print('Ready for M matrix and forcing matrix')
MM_full = kane.mass_matrix_full.subs(kdd).subs(para_dict).subs(in_dict)
F_full = kane.forcing_full.subs(kdd).subs(para_dict).subs(in_dict)


right_side = MM_full.inv()*F_full

######################
# wirting into a file
######################
q_u_dot = array([q1d, q2d, q3d, q4d, u1d, u2d, u3d, u4d, u5d, u6d])

#path
path1 = '/home/stefenyin'
path2 = '/home/stefenstudy'

if os.path.exists(path1):
    path = path1
else:
    path = path2

path = path + '/bicycle/bi_equations_writing/nonlinear_constraints_model_verification.txt'

try:
    f = open(path,'w')
    f.write('')
    f.close()
    del f
    ans = None
    while ans != 'y' and ans != 'n':
        ans = raw_input("%s exists already. Are you sure you want"\
            " to write it? (y or n)\n" % path)
    if ans == 'y':
        f = open(path, 'w')
except IOError:
    f = open(path, 'w')

for i in range(size(q_u_dot)):
    f.write(str(q_u_dot[i]) + '=' + str(right_side[i]) + '\n\n')
