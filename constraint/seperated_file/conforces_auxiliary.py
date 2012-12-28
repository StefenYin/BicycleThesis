###########################################################
# Contact forces in constraint model with axiliary methods
# writing into a file
###########################################################

#==============================================================================
print ('building the model')

import sympy as sym
from numpy import *
import sympy.physics.mechanics as mec

import os

import textwrap

from sympy import signsimp, factor_terms

mec.Vector.simp = False
mec.mechanics_printing()

#-----------------
# Reference Frames


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
        

#-----------------------------------
# Generalized Coordinates and Speeds

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

# Auxiliary generalized speeds for contact point
ua1, ua2, ua3= mec.dynamicsymbols ('ua1 ua2 ua3') 
#for rear wheel contact point, with respect to inertial frame, N1 and N2

ua4, ua5, ua6 = mec.dynamicsymbols ('ua4 ua5 ua6') #for front wheel contact point

ua_zeros_dict = {ua1:0., ua2:0., ua4:0., ua5:0.}


#--------------------------------
# Orientation of Reference Frames

# bicycle frame yaw
A.orient(N, 'Axis', [q1, N['3']])
# bicycle frame roll
B.orient(A, 'Axis', [q2, A['1']])
# bicycle frame pitch
C.orient(B, 'Axis', [q3, B['2']])
# fork/handlebar steer
E.orient(C, 'Axis', [q4, C['3']])

#----------
# parameters

# geometry
# rF: radius of front wheel
# rR: radius of rear wheel
# d1: the perpendicular distance from the steer axis to the center
#    of the rear wheel (rear offset)
# d2: the distance between wheels along the steer axis
# d3: the perpendicular distance from the steer axis to the center
#    of the front wheel (fork offset)
# l1: the distance in the c1> direction from the center of the rear
#    wheel to the frame center of mass
# l2: the distance in the c3> direction from the center of the rear
#    wheel to the frame center of mass
# l3: the distance in the e1> direction from the front wheel center to
#    the center of mass of the fork
# l4: the distance in the e3> direction from the front wheel center to
#    the center of mass of the fork

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

#-------------------
# forces and torques

# control torques
#T2--roll tor, T4--steer tor, T5--rear wheel tor
T4 = mec.dynamicsymbols('T4')

#road_wheel contact forces
Fx_r, Fy_r, Fz_r, Fx_f, Fy_f, Fz_f= mec.dynamicsymbols('Fx_r Fy_r Fz_r Fx_f Fy_f Fz_f')

#------------------------------------      
# Kinematical Differential Equations


kinematical = [q1d - u1,
               q2d - u2,
               q3d - u3,
               q4d - u4]
 
#-------------------
# Angular Velocities


A.set_ang_vel(N, u1 * N['3'])
B.set_ang_vel(A, u2 * A['1'])
C.set_ang_vel(B, u3 * B['2'])
E.set_ang_vel(C, u4 * C['3'])

D.set_ang_vel(C, u5 * C['2'])
F.set_ang_vel(E, u6 * E['2'])
 

print('ready for special unit vectors; and points and velocities')

#--------------------
#special unit vectors

#pitch back for the unit vector g_3 along front wheel radius;

g_3 =  (mec.express(A['3'], E) - mec.dot(E['2'], A['3'])*E['2']).normalize() 
#another way: g_3 = E['2'].cross(A['3']).cross(E['2']).normalize() 


#roll back for longitudinal and lateral unit vector of front wheel

long_v = mec.cross (E['2'], A['3']).normalize()
lateral_v = mec.cross (A['3'], long_v).normalize() 


#---------------------
#points and velocities

###################################################rear wheel contact point, dn
dn = mec.Point('dn')
dn.set_vel(N, ua1 * A['1'] + ua2 * A['2'] + ua3 * A['3'])
#dn.set_vel(N, ua1 * N['1'] + ua2 * N['2'] + ua3 * N['3'])


# rear wheel center, do
do = dn.locatenew('do', -rR * B['3']) 
#do = dn.locatenew('do', -rtr * A['3'] - rR * B['3']) 
# rear wheel center with rear tire radius rtr

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


##################################################Front wheel contact point, fn
fn = mec.Point('fn')
fn.set_vel(N, ua4 * long_v + ua5 * lateral_v + ua6 * A['3'])
#fn.set_vel(N, ua4 * N['1'] + ua5 * N['2'] + ua6 * N['3'])


#front wheel center
fo = fn.locatenew('fo', -rF * g_3)
#fo = fn.locatenew('fo', -ftr * A['3'] - rF * g_3) 
# rear wheel center with rear tire radius rtr.

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
#-------------------
# Motion Constraints
holonomic = [fn.pos_from(dn).dot(A['3'])]

nonholonomic = [(v_SAF_1-v_SAF_2).dot(uv) for uv in E]

#--------
# Inertia

Ic = mec.inertia(C, ic11, ic22, ic33, 0.0, 0.0, ic31)
Id = mec.inertia(C, id11, id22, id11, 0.0, 0.0, 0.0) #rear wheel
Ie = mec.inertia(E, ie11, ie22, ie33, 0.0, 0.0, ie31)
If = mec.inertia(E, if11, if22, if11, 0.0, 0.0, 0.0) #front wheel


#-------------
# Rigid Bodies

rearFrame_inertia = (Ic, co)
rearFrame=mec.RigidBody('rearFrame',co,C,mc,rearFrame_inertia)

rearWheel_inertia = (Id, do)
rearWheel=mec.RigidBody('rearWheel',do,D,md,rearWheel_inertia)

frontFrame_inertia = (Ie, eo)
frontFrame=mec.RigidBody('frontFrame',eo,E,me,frontFrame_inertia)

frontWheel_inertia = (If, fo)
frontWheel=mec.RigidBody('frontWheel',fo,F,mf,frontWheel_inertia)

bodyList = [rearFrame, rearWheel, frontFrame, frontWheel]


#--------------------------
# Generalized Active Forces

# gravity
Fco = (co, mc * g * A['3'])
Fdo = (do, md * g * A['3'])
Feo = (eo, me * g * A['3'])
Ffo = (fo, mf * g * A['3'])

#road_wheel contact forces
F_r = (dn, Fx_r * A['1'] + Fy_r * A['2'] + Fz_r * A['3'])
F_f = (fn, Fx_f * long_v + Fy_f * lateral_v + Fz_f * A['3'])

#F_r = (dn, Fx_r * N['1'] + Fy_r * N['2'] + Fz_r * N['3'])
#F_f = (fn, Fx_f * N['1'] + Fy_f * N['2'] + Fz_f * N['3'])

#input torques
Tc = (C, - T4 * C['3']) #back steer torque to rear frame

Te = (E, T4 * C['3']) #steer torque to front frame


forceList = [Fco, Fdo, Feo, Ffo, F_r, F_f, Tc, Te]


#==============================================================================
print ('Kanes method')


kane = mec.KanesMethod(N, q_ind=[q1, q2, q4], u_ind=[u2, u4, u5], 
    kd_eqs=kinematical, q_dependent=[q3], configuration_constraints=holonomic, 
    u_dependent=[u1, u3, u6], velocity_constraints=nonholonomic, 
    u_auxiliary=[ua1, ua2, ua3, ua4, ua5, ua6])
#reminder: u1--yaw rate, u2--roll rate, u3--pitch rate, u4--steer rate, 
#u5--rear wheel ang. vel., u6--front wheel ang. vel.

(fr, frstar)= kane.kanes_equations(forceList, bodyList)
kdd = kane.kindiffdict()

#path
path1 = '/home/stefenyin'
path2 = '/home/stefenstudy'

if os.path.exists(path1):
    path = path1
else:
    path = path2

con_forces_noncontri = kane.auxiliary_eqs.applyfunc(
        lambda w: factor_terms(signsimp(w))).subs(kdd)

CF = con_forces_noncontri
"""
#==============================================================================
print('building a file and writing the equations into it')

#-----------------
print('building')
path = path + '/bicycle/bi_equations_writing/con_force_nonslip_non_para_input.txt'

try:
    f = open(path,'w')
    f.write('')
    f.close()
    del f
    ans = None
    while ans != 'y' and ans != 'n':
        ans = raw_input("%s exists already. Are you sure you want"\
            " to write contact force equations into"\
            " it? (y or n)\n" % path)
    if ans == 'y':
        f = open(path, 'w')
except IOError:
    f = open(path, 'w')

#----------------
print('writing')
f.write('NOTE:\n' + 
    '1. T4 is an inner torque;\n' +
    '2. Pull force may exist in the experiment, but here in '\
    'the model did not give the variable;\n' +
    '3, No parameters are input, and the contact forces are expressed in '\
    'N inertial frame for each wheel;\n'
    '4. The q3 pitch angle is lam for this model (steer tilt angle).\n\n')


ud_str = ['Derivative(u1(t), t)', 
    'Derivative(u2(t), t)',
    'Derivative(u3(t), t)',
    'Derivative(u4(t), t)',
    'Derivative(u5(t), t)',
    'Derivative(u6(t), t)']

qu_str = ['q1(t)', 'q2(t)', 'q3(t)', 'q4(t)',
    'u1(t)', 'u2(t)', 'u3(t)', 'u4(t)', 'u5(t)', 'u6(t)'] 

CF_str = ['Fx_r(t)', 'Fy_r(t)', 'Fz_r(t)', 'Fx_f(t)', 'Fy_f(t)', 'Fz_f(t)', 'T4(t)']

ud_symbols = ['u1d','u2d','u3d','u4d','u5d','u6d']

qu_symbols = ['q1','q2','q3','q4',
          'u1','u2','u3','u4','u5','u6']

CF_symbols = ['Fx_r', 'Fy_r',  'Fz_r', 'Fx_f', 'Fy_f',  'Fz_f', 'T4']

para_symbols = ['rF','rR', 
        'd1','d2','d3', 
        'l1','l2','l3', 'l4',
        'g', 
        'mc','md','me','mf', 
        'ic11','ic22','ic33','ic31', 
        'id11','id22', 
        'ie11','ie22','ie33','ie31', 
        'if11','if22'
        'v']

wrapper = textwrap.TextWrapper(width=73, initial_indent='   ',
                                subsequent_indent='        ')

number = len(CF)

for num in range(number):
    exp = CF[num].evalf(n=3)
    CF[num] = exp

    exp_str = str(exp)

    title = 'Expression for: '

    para = 'Parameters: '

    signals = 'Dynamic symbols: '

    force = 'Contact forces symbols: ' 

    for a,b in zip(ud_symbols, ud_str):
        if b in exp_str:
            signals += a + ', '

        exp_str = exp_str.replace(b, a)
   
    for c,d in zip(qu_symbols, qu_str):
        if d in exp_str:
            signals += c + ', '

        exp_str = exp_str.replace(d, c)
        
    for e,h in zip(CF_symbols, CF_str):
        if h in exp_str:
            title += e + ', '
            force += e + ', '
        exp_str = exp_str.replace(h, e)

    for i in para_symbols:
        if i in exp_str:
            para += i + ', '
    
    if 't' in exp_str:
        note = 'It still has some variables in terms of <t>'
    else:
        note = 'No variables about <t> at all'

    f.write (title + '\n')
    f.write (para + '\n')
    f.write (signals + '\n')
    f.write (force + '\n')
    f.write (note + '.\n\n')

    exp_str = wrapper.wrap(exp_str)
    for k in range(len(exp_str)):
        f.write (exp_str[k] + ' \\' + '\n')

    f.write ('\n\n')

f.close()
"""
