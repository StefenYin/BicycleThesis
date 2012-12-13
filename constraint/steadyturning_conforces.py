#############################################################
# steady turning simulation when given lean and steer angles.
# for contact forces
#############################################################

#################
# nonlinear model
#################

print ("Nonlinear bicycle whipple model.")

import sympy as sym
from numpy import *
from numpy import sin, cos, sqrt, arctan
import sympy.physics.mechanics as mec

import os

from dtk import bicycle as bi
mec.Vector.simp = False
mec.mechanics_printing()

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


# bicycle frame yaw
A.orient(N, 'Axis', [q1, N['3']])
# bicycle frame roll
B.orient(A, 'Axis', [q2, A['1']])
# bicycle frame pitch
C.orient(B, 'Axis', [q3, B['2']])
# fork/handlebar steer
E.orient(C, 'Axis', [q4, C['3']])

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

# control torques
#T2--roll tor, T4--steer tor, T5--rear wheel tor
T4 = mec.dynamicsymbols('T4')

kinematical = [q1d - u1,
               q2d - u2,
               q3d - u3,
               q4d - u4]

A.set_ang_vel(N, u1 * N['3'])
B.set_ang_vel(A, u2 * A['1'])
C.set_ang_vel(B, u3 * B['2'])
E.set_ang_vel(C, u4 * C['3'])

D.set_ang_vel(C, u5 * C['2'])
F.set_ang_vel(E, u6 * E['2'])


#pitch back for the unit vector g_3 along front wheel radius;

g_3 =  (mec.express(A['3'], E) - mec.dot(E['2'], A['3'])*E['2']).normalize() 
#another way: g_3 = E['2'].cross(A['3']).cross(E['2']).normalize() 


#roll back for longitudinal and lateral unit vector of front wheel

long_v = mec.cross (E['2'], A['3']).normalize()
lateral_v = mec.cross (A['3'], long_v).normalize() 

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

holonomic = [fn.pos_from(dn).dot(A['3'])]

nonholonomic = [(v_SAF_1-v_SAF_2).dot(uv) for uv in E]

Ic = mec.inertia(C, ic11, ic22, ic33, 0.0, 0.0, ic31)
Id = mec.inertia(C, id11, id22, id11, 0.0, 0.0, 0.0) #rear wheel
Ie = mec.inertia(E, ie11, ie22, ie33, 0.0, 0.0, ie31)
If = mec.inertia(E, if11, if22, if11, 0.0, 0.0, 0.0) #front wheel

rearFrame_inertia = (Ic, co)
rearFrame=mec.RigidBody('rearFrame',co,C,mc,rearFrame_inertia)

rearWheel_inertia = (Id, do)
rearWheel=mec.RigidBody('rearWheel',do,D,md,rearWheel_inertia)

frontFrame_inertia = (Ie, eo)
frontFrame=mec.RigidBody('frontFrame',eo,E,me,frontFrame_inertia)

frontWheel_inertia = (If, fo)
frontWheel=mec.RigidBody('frontWheel',fo,F,mf,frontWheel_inertia)

bodyList = [rearFrame, rearWheel, frontFrame, frontWheel]

# gravity
Fco = (co, mc * g * A['3'])
Fdo = (do, md * g * A['3'])
Feo = (eo, me * g * A['3'])
Ffo = (fo, mf * g * A['3'])

# input torques
Tc = (C, - T4 * C['3']) #back steer torque to rear frame

Te = (E, T4 * C['3']) #steer torque to front frame

forceList = [Fco, Fdo, Feo, Ffo, Tc, Te]

kane = mec.KanesMethod(N, q_ind=[q1, q2, q4], u_ind=[u2, u4, u5], kd_eqs=kinematical, 
                       q_dependent=[q3], configuration_constraints=holonomic, 
                       u_dependent=[u1, u3, u6], velocity_constraints=nonholonomic)
#reminder: u1--yaw rate, u2--roll rate, u3--pitch rate, u4--steer rate, u5--rear wheel ang. vel., u6--front wheel ang. vel.

#kane.kindiffeq(kinematical)
(fr, frstar)= kane.kanes_equations(forceList, bodyList)
kdd = kane.kindiffdict()

print('Ready for nonlinear model')
MM_full = kane.mass_matrix_full.subs(kdd)
F_full = kane.forcing_full.subs(kdd)


#################################
# steady turning equation building
###################################

#==============================================================================
print ("pitch angle function from holonomic equation.")

from scipy.optimize import newton

def pitch_from_roll_and_steer(q2, q4, rF, rR, d1, d2, d3, guess=None):
    def pitch_constraint(q3, q2, q4, rF, rR, d1, d2, d3):
        zero = -d1*sin(q3)*cos(q2) - rR*cos(q2) + \
        (d2 + rF*cos(q2)*cos(q3)/sqrt((sin(q2)*sin(q4) - 
        sin(q3)*cos(q2)*cos(q4))**2 + 
        cos(q2)**2*cos(q3)**2))*cos(q2)*cos(q3) + \
        (d3 + rF*(sin(q2)*sin(q4) - sin(q3)*cos(q2)*
        cos(q4))/sqrt((sin(q2)*sin(q4) - 
        sin(q3)*cos(q2)*cos(q4))**2 + cos(q2)**2*
        cos(q3)**2))*(sin(q2)*sin(q4) - sin(q3)*cos(q2)*cos(q4))
        
        return zero

    if guess is None:
        # guess based on steer and roll being both zero
        guess = bi.lambda_from_abc(rF, rR, d1, d3, d2)

    args = (q2, q4, rF, rR, d1, d2, d3)

    q3 = newton(pitch_constraint, guess, args=args)

    return q3

#==============================================================================
print ('parameters and states values')

#benchmark parameters from bicycle module
bp = bi.benchmark_parameters()

mp = bi.benchmark_to_moore(bp)

para_dict={l1: mp['l1'], l2: mp['l2'], l3: mp['l3'], l4:mp['l4'], 
        mc: mp['mc'], md:mp['md'], me:mp['me'], mf:mp['mf'], rR:mp['rr'], rF:mp['rf'], 
        g:mp['g'], d1: mp['d1'], d2: mp['d2'], d3: mp['d3'],
        id11: mp['id11'], id22: mp['id22'], if11: mp['if11'], if22: mp['if22'],
        ic11: mp['ic11'], ic22: mp['ic22'], ic33: mp['ic33'], ic31: mp['ic31'],
        ie11: mp['ie11'], ie22: mp['ie22'], ie33: mp['ie33'], ie31: mp['ie31']}

#steady turning configuration
#left side states: q1d, q2d, q4d, q3d, u2d, u4d, u5d, u1d, u3d, u6d. 
#Only yaw rate left.
q_u_dot = array([u1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#assign the roll rate, pitch rate, and steer rate to be zeros
u_dict = {u2: 0.0, u3: 0.0, u4: 0.0}


#leansteer and pitch
lean = pi/8
steer = pi/4

#If you know the limits, please use the try-except-else for pitch angle.
pitch = bi.pitch_from_roll_and_steer(lean, steer, mp['rf'], mp['rr'], mp['d1'], 
        mp['d2'], mp['d3'])
# pitch = pitch_from_roll_and_steer(lean, steer, mp['rf'], mp['rr'], mp['d1'], 
#        mp['d2'], mp['d3'])

q_dict = {q2: lean, q3:pitch, q4:steer}

#==============================================================================
print ('three dynamics equations')

print ("a*(u1*u1) + b*(u1*u5) + c*(u1*u6) + d*T4 + e = 0\n")

F_full_1 = F_full.subs(u_dict).subs(para_dict).subs(q_dict).expand()

dynamic_equ = F_full_1[4:7]
#4:7 means the dynamic equations, 
# we need the nonholonomic questions for solve the values.

num_dict = {u1: 0.0, u5: 0.0, u6:0.0, T4: 0.0}
dynamic_equ_coeff = [[value.coeff(u1*u1), value.coeff(u1*u5), value.coeff(u1*u6), 
            value.coeff(T4)] for value in dynamic_equ]

dynamic_equ_indi = [-value.subs(num_dict) for value in dynamic_equ]

#==============================================================================
print ('nonholonomic equations')

nonholonomic_2 = [value.subs(para_dict).subs(q_dict) for value in nonholonomic]

nonholonomic_coeff_inde = [[value.coeff(u2), value.coeff(u4), 
                value.coeff(u5)] for value in nonholonomic_2]

nonholonomic_coeff_de = [[value.coeff(u1), value.coeff(u3), 
                value.coeff(u6)] for value in nonholonomic_2]


inde_states = matrix([[u2], [u4], [u5]])
nonho_coeff_inde_ma = asmatrix(nonholonomic_coeff_inde)
nonho_coeff_de_ma = asmatrix(nonholonomic_coeff_de)

left_side_ma = (nonho_coeff_de_ma.I) * (- nonho_coeff_inde_ma * inde_states)
left_side_arr = asarray(left_side_ma)

left_side_subs_list = [value[0].subs(u_dict) for value in left_side_arr]
u1_equ = left_side_subs_list[0]
u6_equ = left_side_subs_list[2]

print u1_equ, u6_equ

#==============================================================================
print ('substitute the nonholonomic equations into dynamic euqations.')
nonho_dict = dict(zip([u1, u6], [u1_equ, u6_equ])) #expression for u1 and u6
dynamic_equ_subs = [value.subs(nonho_dict).expand() for value in dynamic_equ] 
#two unknowns u5 and T4 in three equations

#calculate the u5, T4 and resulting u1 and u6
u5_value = sym.solve(dynamic_equ_subs[0], u5)[0]  #choose the negative value
nonho_dict.update({u5: u5_value})

T4_value = sym.solve(dynamic_equ_subs[1], T4)[0].subs(nonho_dict)

u6_value = nonho_dict[u6].subs(nonho_dict)
u1_value = nonho_dict[u1].subs(nonho_dict)

nonho_dict_value = dict(nonho_dict) #just cope it
nonho_dict_value[u6] =u6_value
nonho_dict_value[u1] = u1_value

print nonho_dict_value  #u1, u5 and u6
print T4_value
