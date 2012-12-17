####################################
# linearized benchmark model
# check with benchmark eigenvalues
####################################

import sympy as sym
from numpy import *
import sympy.physics.mechanics as mec

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
#do = dn.locatenew('do', -rtr * A['3'] - rR * B['3']) # rear wheel center 
#with rear tire radius rtr

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
#fo = fn.locatenew('fo', -ftr * A['3'] - rF * g_3) # rear wheel center with 
#rear tire radius rtr.

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

kane = mec.KanesMethod(N, q_ind=[q1, q2, q4], u_ind=[u2, u4, u5], 
    kd_eqs=kinematical, q_dependent=[q3], configuration_constraints=holonomic, 
    u_dependent=[u1, u3, u6], velocity_constraints=nonholonomic)
#reminder: u1--yaw rate, u2--roll rate, u3--pitch rate, u4--steer rate, 
#u5--rear wheel ang. vel., u6--front wheel ang. vel.

#kane.kindiffeq(kinematical)
(fr, frstar)= kane.kanes_equations(forceList, bodyList)
kdd = kane.kindiffdict()


############################
#substituting numerical values
#############################
#parameters from Meijaard, paper 2007
w = 1.02                     # wheelbase         [m]
c = 0.08                     # trail             [m]
lam = pi/10                    # steer axis tilt   [rad]
S_g = 9.81                     # gravity           [N/kg]

#rear wheel
S_rR     =   0.30                     # rear wheel radius [m]
mR     =   2.00                     # rear wheel mass   [kg] 
IRxx   =   0.0603                   #                   [kg*m^2]
IRyy   =   0.12                     #                   [kg*m^2]

#rear frame
xB     =   0.30                     #                   [m]
zB     = - 0.90                     #                   [m]
mB     =  85.00                     # frame mass        [kg]
IBxx   =   9.20                     #                   [kg*m^2]
IBxz   =   2.40                     #                   [kg*m^2]
IByy   =  11.00                     #                   [kg*m^2]
IBzz   =   2.80                     #                   [kg*m^2]

#front frame
xH     =   0.90                     #                   [m]
zH     = - 0.70                     #                   [m]
mH     =   4.00                     #                   [kg]
IHxx   =   0.05892                  #                   [kg*m^2]
IHxz   = - 0.00756                  #                   [kg*m^2]
IHyy   =   0.06                     #                   [kg*m^2]
IHzz   =   0.00708                  #                   [kg*m^2]

#front wheel
S_rF     =   0.35                     #                   [m]
mF     =   3.00                     #                   [kg]
IFxx   =   0.1405                   #                   [kg*m^2]
IFyy   =   0.28                     #                   [kg*m^2]

#rotation of the fork
IH     = mat([[IHxx,  0,  IHxz],
          [0,    IHyy, 0],   
          [IHxz, 0,    IHzz]])

cf     = mat([[cos(lam), 0, -sin(lam)],
          [0, 1,  0],          
          [sin(lam), 0,  cos(lam)]])
          
IHrot  = cf*IH*(cf.T)         

#rotation of the front frame
IB     = mat([[IBxx, 0, IBxz],
             [0,  IByy,  0],
             [IBxz,  0,  IBzz]])
             
IBrot  = cf*IB*(cf.T)

#parameters for Imodel transfered from 2007 paper
S_d1   =  cos(lam)*(c+w-S_rR*tan(lam))
S_d2   =  -cos(lam)*(S_rF-S_rR-w*tan(lam))
S_d3   =  -cos(lam)*(c-S_rF*tan(lam))

S_l1   =  xB*cos(lam)-zB*sin(lam)-S_rR*sin(lam)
S_l2   =  xB*sin(lam)+zB*cos(lam)+S_rR*cos(lam)
S_l4   =  (zH+S_rF)*cos(lam)+(xH-w)*sin(lam)
S_l3   =  -(zH+S_rF)*sin(lam)+(xH-w)*cos(lam)
#S_l3  =  (xH-w-S_l4*sin(lam))/cos(lam) #this one has association with S_l4

#rear wheel. Here, the radius is S_rR
S_id11                            =  IRxx                   
S_id22                            =  IRyy                    
S_id33                            =  IRxx                    
S_md                              =  mR

#rear frame
S_ic11                            =  IBrot[0,0]                     
S_ic12                            =  IBrot[0,1]                     
S_ic22                            =  IBrot[1,1]                      
S_ic23                            =  IBrot[1,2]                      
S_ic31                            =  IBrot[2,0]                      
S_ic33                            =  IBrot[2,2]                      
S_mc                              =  mB  

#front frame
S_ie11                            =  IHrot[0,0]              
S_ie12                            =  IHrot[0,1]              
S_ie22                            =  IHrot[1,1]              
S_ie23                            =  IHrot[1,2]              
S_ie31                            =  IHrot[2,0]                             
S_ie33                            =  IHrot[2,2]                             
S_me                              =  mH

#front wheel. Here the radius is S_rF
S_if11                            =  IFxx                                   
S_if22                            =  IFyy                                   
S_if33                            =  IFxx                                  
S_mf                              =  mF
                                                                                
# speed v and build the dictionary 
v = sym.Symbol('v')

para_dict={d1: S_d1, d2: S_d2, d3: S_d3, id11: S_id11, id22: S_id22, 
    ic11: S_ic11, ic22: S_ic22, ic31: S_ic31, ic33: S_ic33, ie11: S_ie11, 
    ie22: S_ie22, ie31: S_ie31, ie33: S_ie33, if11: S_if11, if22: S_if22, 
    l1: S_l1, l2: S_l2, l3: S_l3, l4:S_l4, md:S_md, mc:S_mc, me:S_me, mf:S_mf, 
    rR:S_rR, rF:S_rF, g:S_g} #for nonlinear substitutioin

qu_dict = { q1: 0., q2: 0., q3:lam, q4:0., u1:0., u2:0., u3:0., u4:0., 
        u5: -v/S_rR, u6: -v/S_rF}


print('Ready for M matrix')
MM_full = kane.mass_matrix_full
MM_full_1 = MM_full.subs(para_dict).subs(qu_dict)

print('Ready for Linearization of the Forcing term and B matrix')
forcing_lin_A= kane.linearize()[0].subs(kdd)
forcing_lin_A_1 = forcing_lin_A.subs(para_dict).subs(qu_dict)

#forcing_lin_B = kane.linearize()[1].subs(kdd)
#forcing_lin_B_1 = forcing_lin_B.subs(para_dict).subs(qu_dict)

print('ready for evalf() call of MM_full_2 and forcing_lin')
MM_full_2 = MM_full_1.evalf()
forcing_lin_A_2 = forcing_lin_A_1.evalf()
#forcing_lin_B_2 = forcing_lin_B_1.evalf()

print('ready for A matrix')
Amat = MM_full_2.inv() * forcing_lin_A_2
#Bmat = MM_full_2.inv() * forcing_lin_B_2

Am = Amat.extract([1,2,4,5],[1,2,3,4])
