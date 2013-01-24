"""
In this module, I am going to derive the model mannually with Python codes
using Kane's method in order to calculate the contact forces.
Good Luck!
"""

import sympy.physics.mechanics as mec
import sympy as sym
import numpy as np
import os
import textwrap

from sympy import (sin, cos, sqrt)
from time import time
from writing import writing_dyequ

np.set_printoptions(precision=17)

mec.Vector.simp = False

# Generalized Coordinates and Speeds:
# q1,u1: frame yaw
# q2,u2: frame roll
# q3,u3: frame pitch 
# q4,u4: steering rotation
# u5: rear wheel ang. vel.
# u6: front wheel ang. vel.
# Axiliary speeds at contact points:
# Rear wheel: ua1, ua2
# Front wheel: ua4, ua5
# q = [q_in + q_de]
# u = [u_in + u_de]
q1, q2, q3, q4 = mec.dynamicsymbols('q1 q2 q3 q4')
q1d, q2d, q3d, q4d = mec.dynamicsymbols('q1 q2 q3 q4', 1)
u1, u2, u3, u4 = mec.dynamicsymbols('u1 u2 u3 u4')
u5, u6 = mec.dynamicsymbols('u5 u6')
u1d, u2d, u3d, u4d = mec.dynamicsymbols('u1 u2 u3 u4', 1)
u5d, u6d = mec.dynamicsymbols('u5 u6', 1)
ua1, ua2, ua3 = mec.dynamicsymbols ('ua1 ua2 ua3')
ua4, ua5, ua6 = mec.dynamicsymbols ('ua4 ua5 ua6')

#
ua1d, ua2d, ua3d = mec.dynamicsymbols ('ua1 ua2 ua3',1)
ua4d, ua5d, ua6d = mec.dynamicsymbols ('ua4 ua5 ua6',1)

q = [q1, q2, q4, q3]
qd = [q1d, q2d, q4d, q3d]
u = [u2, u4, u5, u1, u3, u6]
ud = [u2d, u4d, u5d, u1d, u3d, u6d]
ua = [ua1, ua2, ua3, ua4, ua5, ua6]
#
uad = [ua1d, ua2d, ua3d, ua4d, ua5d, ua6d]

q_zero = {qi: 0 for qi in q}
qd_zero = {qdi: 0 for qdi in qd}
u_zero = {ui: 0 for ui in u}
ud_zero = {udi: 0 for udi in ud}
ua_zero = {uai: 0 for uai in ua}
#
uad_zero = {uadi: 0 for uadi in uad}

# Active Forces:
# T4: steer torque.
# Fx_r, Fy_r, Fz_r: rear wheel contact forces.
# Fx_f, Fy_f, Fz_f: front wheel contact forces.
# Fco, Fdo, Feo, Ffo: gravity of four rigid bodies of bicycle.
T4 = mec.dynamicsymbols('T4')
Fx_r, Fy_r, Fz_r, Fx_f, Fy_f, Fz_f = mec.dynamicsymbols('Fx_r Fy_r Fz_r Fx_f Fy_f Fz_f')

auxforces = [Fx_r, Fy_r, Fz_r, Fx_f, Fy_f, Fz_f]

# Parameters:
# Geometry:
# rf: radius of front wheel
# rr: radius of rear wheel
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
# Mass:
# mc, md, me, mf: mass of rearframe, rearwheel, frontframe, frontwheel
# Inertia:
# ic11, ic22, ic33, ic31: rear frame
# id11, id22: rear wheel
# ie11, ie22, ie33, ie31: front frame
# if11, if22: front wheel
rf, rr = sym.symbols('rf rr')
d1, d2, d3 = sym.symbols('d1 d2 d3')
l1, l2, l3, l4 = sym.symbols('l1 l2 l3 l4')

g = sym.symbols('g')
t = sym.symbols('t')

mc, md, me, mf = sym.symbols('mc md me mf')

ic11, ic22, ic33, ic31 = sym.symbols('ic11 ic22 ic33 ic31')
id11, id22 = sym.symbols('id11 id22')
ie11, ie22, ie33, ie31 = sym.symbols('ie11 ie22 ie33 ie31')
if11, if22 = sym.symbols('if11 if22')

# Reference Frames:
# Newtonian Frame: N
# Yaw Frame: A
# Roll Frame: B
# Pitch & Bicycle Frame: C
# Steer & Fork/Handlebar Frame: E
# Rear Wheel Frame: D
# Front Wheel Frame: F
# Orientation of Reference Frames:
# bicycle frame yaw: N->A
# bicycle frame roll: A->B
# pitch to rear frame: B->C
# fork/handlebar steer: C->E
N = mec.ReferenceFrame('N', indices=('1', '2', '3'))
A = mec.ReferenceFrame('A', indices=('1', '2', '3'))
B = mec.ReferenceFrame('B', indices=('1', '2', '3'))
C = mec.ReferenceFrame('C', indices=('1', '2', '3'))
E = mec.ReferenceFrame('E', indices=('1', '2', '3'))
D = mec.ReferenceFrame('D', indices=('1', '2', '3'))
F = mec.ReferenceFrame('F', indices=('1', '2', '3'))

A.orient(N, 'Axis', [q1, N['3']])
B.orient(A, 'Axis', [q2, A['1']])
C.orient(B, 'Axis', [q3, B['2']])
E.orient(C, 'Axis', [q4, C['3']])

# Angular Velocities
# Define the generalized speeds:
A.set_ang_vel(N, u1 * N['3'])
B.set_ang_vel(A, u2 * A['1'])
C.set_ang_vel(B, u3 * B['2'])
E.set_ang_vel(C, u4 * C['3'])
D.set_ang_vel(C, u5 * C['2'])
F.set_ang_vel(E, u6 * E['2'])

# Kinematical Differential Equations:
kindiffs = sym.Matrix([q1d - u1,
                       q2d - u2,
                       q3d - u3,
                       q4d - u4])
qd_kd = sym.solve(kindiffs, qd)

# Values of generalized speeds during a steady turning
steady_conditions = sym.solve(kindiffs.subs({q2d:0, q3d:0, q4d:0}), [u1, u2, u3, u4])
steady_conditions.update({q2d:0, q3d:0, q4d:0})

# Special unit vectors:
# g_3: direction along front wheel radius.
# long_v: longitudinal direction of front wheel.
# lateral_v: lateral direction of front wheel.
g_3 =  (mec.express(A['3'], E) - mec.dot(E['2'], A['3'])*E['2']).normalize() 
# OR g_3 = E['2'].cross(A['3']).cross(E['2']).normalize() 
long_v = mec.cross (E['2'], A['3']).normalize()
lateral_v = mec.cross (A['3'], long_v)

# Points and velocities:
# dn: rear wheel contact point.
# do: rear wheel center.
# rtr: rear tire radius.
# fn: front wheel contact point.
# fo: front wheel center.
# ftr: front tire radius
# co: rear frame center.
# eo: front frame center.
# ce: steer axis point.
# SAF: steer axis foot.
## Rear
dn = mec.Point('dn')
dn.set_vel(N, ua1 * A['1'] + ua2 * A['2'] + ua3 * A['3']) 

do = dn.locatenew('do', (-rr * B['3']).express(C))
do.v2pt_theory(dn, N, D)
do.set_acc(N, do.vel(N).subs(ua_zero).diff(t, C).subs(qd_kd) + 
           mec.cross(C.ang_vel_in(N), do.vel(N).subs(ua_zero).subs(qd_kd)))

co = do.locatenew('co', l1 * C['1'] + l2 * C['3'])
co.v2pt_theory(do, N, C)
#co.a2pt_theory(do, N, C)
co.set_acc(N, co.vel(N).subs(ua_zero).diff(t, C).subs(qd_kd) + 
           mec.cross(C.ang_vel_in(N), co.vel(N).subs(ua_zero).subs(qd_kd)))

## Front
fn = mec.Point('fn')
fn.set_vel(N, ua4 * long_v + ua5 * lateral_v + ua6 * A['3'])

fo = fn.locatenew('fo', (-rf * g_3).express(E))
fo.v2pt_theory(fn, N, F)
fo.set_acc(N, fo.vel(N).subs(ua_zero).diff(t, E).subs(qd_kd) + 
           mec.cross(E.ang_vel_in(N), fo.vel(N).subs(ua_zero).subs(qd_kd)))

eo = fo.locatenew('eo', l3 * E['1'] + l4 * E['3'])
eo.v2pt_theory(fo, N, E)
#eo.a2pt_theory(fo, N, E)
eo.set_acc(N, eo.vel(N).subs(ua_zero).diff(t, E).subs(qd_kd) + 
           mec.cross(E.ang_vel_in(N), eo.vel(N).subs(ua_zero).subs(qd_kd)))

# Holo and nonholo Constraints
# f_v: B1(q,t) * [u] + B2(t) [ua] + B3= 0
# f_a: B4(q,t) * [du/dt] + B5(q, dq/dt, u, t) + B6(t) = 0
# where, B3 and B6 are zeros; B1 == B4; zeros coefficients of dq/dt.
SAF = do.locatenew('SAF', d1 * C['1'] + d2 * C['3'])
SAF.set_pos(fo, -d3 * E['1'])

v_SAF_1 = do.vel(N) + mec.cross(C.ang_vel_in(N), SAF.pos_from(do))
a_SAF_1 = do.acc(N) + mec.cross(C.ang_vel_in(N), mec.cross(C.ang_vel_in(N), SAF.pos_from(do))) + \
          mec.cross(C.ang_acc_in(N), SAF.pos_from(do))

v_SAF_2 = fo.vel(N) + mec.cross(E.ang_vel_in(N), SAF.pos_from(fo))
a_SAF_2 = fo.acc(N) + mec.cross(E.ang_vel_in(N), mec.cross(E.ang_vel_in(N), SAF.pos_from(fo))) + \
          mec.cross(E.ang_acc_in(N), SAF.pos_from(fo))

f_c = sym.Matrix([fn.pos_from(dn).dot(A['3'])])
f_v = sym.Matrix([(v_SAF_1 - v_SAF_2).dot(uv) for uv in E])
f_a = sym.Matrix([(a_SAF_1 - a_SAF_2).dot(uv) for uv in E])

"""
# This way is gonna to used for linearization if possible.
# Built B1, B2, B3 and assert
B3 = f_v.subs(u_zero).subs(ua_zero)
B1 = sym.zeros(3, len(u))
B2 = sym.zeros(3, len(ua))
B4 = sym.zeros(3, len(u))
B5 = sym.zeros(3, 1)

for i, (f_vi, f_ai) in enumerate(zip(f_v, f_a)):
    for j, (ui, udi) in enumerate(zip(u, ud)):
        B1[i, j] = f_vi.diff(ui)
        B4[i, j] = f_ai.diff(udi)
        assert B1[i, j] == B4[i, j]
    for k, uai in enumerate(ua):
        B2[i, k] = f_vi.diff(uai)
    B5[i] = f_ai.subs(ud_zero)
    assert B5[i].subs(u_zero).subs(q_zero) == 0
    for qdi in qd:
        assert B5[i].diff(qdi) == 0
"""

# Acceleration of constraints checks.
B4 = sym.zeros(3, len(u))
B5 = sym.zeros(3, 1)
for i, f_ai in enumerate(f_a):
    for j, udi in enumerate(ud):
        B4[i, j] = f_ai.diff(udi)
    B5[i] = f_ai.subs(ud_zero)
    assert B5[i].subs(u_zero).subs(q_zero).subs(ua_zero) == 0
    for qdi in qd:
        assert B5[i].diff(qdi).doit() == 0


t0 = time()
# Constraint coefficient matrix: M_v * [u; ua] = 0
# assert B1.row_join(B2) == M_v
M_v = sym.zeros((3, 12))
for i in range(3):
    for j, ui in enumerate(u + ua):
        M_v[i, j] = f_v[i].diff(ui)
# Form above constraint expression into u_de = A_rs * [u_in; ua]
M_v_i = M_v[:, :3]  # take u2, u4, u5
M_v_d = M_v[:, 3:6] # take u1, u3, u6
M_v_aux = M_v[:, 6:] # rest 6 colums for auxiliary
M_v_i_aux = M_v_i.row_join(M_v_aux)

#A_rs = - M_v_d.inv() * M_v_i_aux # this method takes a long time to calculate.
A_rs = -(M_v_d.adjugate() *  M_v_i_aux) / M_v_d.berkowitz_det()

# Depend speeds expressed by independ speeds since ua == 0
u_dep = A_rs[:, :3] * sym.Matrix(u[:3])
u_dep_dict = {udei : u_depi[0] for udei, u_depi in zip(u[3:], u_dep.tolist())}
t1 = time()
print ("Constraint equations for u_de = A_rs * [u_in; ua], taking {0} s".format(
        t1-t0))


# Partial velocities of points and angular velocities of bodies
partial_v_do = [do.vel(N).diff(ui, N) for ui in u + ua]
partial_v_co = [co.vel(N).diff(ui, N) for ui in u + ua]
partial_v_eo = [eo.vel(N).diff(ui, N) for ui in u + ua]
partial_v_fo = [fo.vel(N).diff(ui, N) for ui in u + ua]
partial_v_dn = [dn.vel(N).diff(ui, N) for ui in u + ua]
partial_v_fn = [fn.vel(N).diff(ui, N) for ui in u + ua]

partial_w_D = [D.ang_vel_in(N).diff(ui, N) for ui in u + ua]
partial_w_C = [C.ang_vel_in(N).diff(ui, N) for ui in u + ua]
partial_w_E = [E.ang_vel_in(N).diff(ui, N) for ui in u + ua]
partial_w_F = [F.ang_vel_in(N).diff(ui, N) for ui in u + ua]


t2 = time()
# Active forces and torques
F_do = md * g * A['3']
F_co = mc * g * A['3']
F_eo = me * g * A['3']
F_fo = mf * g * A['3']
F_dn = Fx_r * A['1'] + Fy_r * A['2'] + Fz_r * A['3']
F_fn = Fx_f * long_v + Fy_f * lateral_v + Fz_f * A['3']
T_E = T4 * E['3']
T_C = -T4 * C['3']

# Generalized active forces
# Simplify some Fr_u equations; before it I tested the simplification to make 
# sure it would take short time. Finally, I simplified [3, -2, -3], primarily
# the contact forces of front contact. So you will find the auxiliary forces on
# each wheel appear on each Fr_u equation of the wheel, respectively and 
# individually.
Fr_u = sym.Matrix([mec.dot(F_do, pv_do) + mec.dot(F_co, pv_co) + 
                   mec.dot(F_eo, pv_eo) + mec.dot(F_fo, pv_fo) + 
                   mec.dot(F_dn, pv_dn) + mec.dot(F_fn, pv_fn) + 
                   mec.dot(T_E, pw_E) + mec.dot(T_C, pw_C)
                   for pv_do, pv_co, pv_eo, pv_fo, pv_dn, pv_fn, pw_E, pw_C in
                   zip(partial_v_do, partial_v_co, partial_v_eo, partial_v_fo,
                   partial_v_dn, partial_v_fn, partial_w_E, partial_w_C)])
Fr_u[3] = Fr_u[3].simplify()
Fr_u[-2] = Fr_u[-2].simplify()
Fr_u[-3] = Fr_u[-3].simplify()
t3 = time()
print ("Generalized active forces, taking {0} s".format(t3-t2))


# Inertia
# Inertia forces and torques
I_D = mec.inertia(C, id11, id22, id11)
I_C = mec.inertia(C, ic11, ic22, ic33, 0., 0., ic31)
I_E = mec.inertia(E, ie11, ie22, ie33, 0., 0., ie31)
I_F = mec.inertia(E, if11, if22, if11)

R_star_do = -md * do.acc(N)
R_star_co = -mc * co.acc(N)
R_star_eo = -me * eo.acc(N)
R_star_fo = -mf * fo.acc(N)
T_star_D = -(mec.dot(I_D, D.ang_acc_in(N)) + mec.cross(D.ang_vel_in(N), 
                                             mec.dot(I_D, D.ang_vel_in(N))))
T_star_C = -(mec.dot(I_C, C.ang_acc_in(N)) + mec.cross(C.ang_vel_in(N), 
                                             mec.dot(I_C, C.ang_vel_in(N))))
T_star_E = -(mec.dot(I_E, E.ang_acc_in(N)) + mec.cross(E.ang_vel_in(N), 
                                             mec.dot(I_E, E.ang_vel_in(N))))
T_star_F = -(mec.dot(I_F, F.ang_acc_in(N)) + mec.cross(F.ang_vel_in(N), 
                                             mec.dot(I_F, F.ang_vel_in(N))))

R_star_do_u1 = sym.Matrix([mec.dot(R_star_do, pv_do)
                         for pv_do in partial_v_do[:6]])
R_star_do_u2 = sym.Matrix([mec.dot(R_star_do, pv_do).simplify()
                         for pv_do in partial_v_do[6:]])
R_star_do_u = R_star_do_u1.col_join(R_star_do_u2)
R_star_co_u1 = sym.Matrix([mec.dot(R_star_co, pv_co)
                         for pv_co in partial_v_co[:6]])
R_star_co_u2 = sym.Matrix([mec.dot(R_star_co, pv_co).simplify()
                         for pv_co in partial_v_co[6:]])
R_star_co_u = R_star_co_u1.col_join(R_star_co_u2)
R_star_eo_u = sym.Matrix([mec.dot(R_star_eo, pv_eo) for pv_eo in partial_v_eo])
R_star_fo_u = sym.Matrix([mec.dot(R_star_fo, pv_fo) for pv_fo in partial_v_fo])

T_star_u = sym.Matrix([mec.dot(T_star_D, pw_D) + mec.dot(T_star_C, pw_C) +
                        mec.dot(T_star_E, pw_E) + mec.dot(T_star_F, pw_F)
                        for pw_D, pw_C, pw_E, pw_F
                        in zip(partial_w_D, partial_w_C, partial_w_E,
                        partial_w_F)])
Fr_star_u = R_star_do_u + R_star_co_u + R_star_eo_u + R_star_fo_u + T_star_u
t4 = time()
print ("Generalized inertia forces, taking {0} s".format(t4-t3))


# With A_rs, Fr_u and Fr_star_u -> Fr_c, Fr_star_c
# dynamic equations and forces equations for contact forces.
# dynamic equations form B7(q, t) * du/dt + B8(T4, u, q, t) = 0
Fr_c = Fr_u[:3, :].col_join(Fr_u[6:, :]) + A_rs.T * Fr_u[3:6, :]
Fr_star_c = Fr_star_u[:3, :].col_join(Fr_star_u[6:, :]) + \
            A_rs.T * Fr_star_u[3:6, :]
dynamic = (Fr_c + Fr_star_c)[:3, :]
forces = (Fr_c + Fr_star_c)[3:, :]

B7 = sym.zeros(3, len(u))
B8 = sym.zeros(3, 1)
for i, dyi in enumerate(dynamic):
    for j, udi in enumerate(ud):
        B7[i, j] = dyi.diff(udi)
    B8[i] = dyi.subs(ud_zero)
    for qdi in qd:
        assert B8[i].diff(qdi) == 0


# Import parameters
from bicycle import (benchmark_parameters, benchmark_to_moore)
from steadyturning import configuration
from numpy import pi

bp = benchmark_parameters()
mp = benchmark_to_moore(bp)
para_dict = {}
for key, value in mp.items():
    para_dict.update(dict(zip([sym.symbols(key)], [value])))


"""
# Validate the nonlinear model
from bicycle import (basu_table_one_input, basu_to_stefen_input, 
                     basu_table_one_output, basu_to_stefen_output)

T4_zero = {T4: 0.}

basu_input = basu_table_one_input()
stefen_input = basu_to_stefen_input(basu_input, mp['rr'], bp['lambda'])
input_states_dict = {}
for key, value in stefen_input.items():
    input_states_dict.update(dict(zip([mec.dynamicsymbols(key)], [value])))

basu_output = basu_table_one_output()
stefen_output = basu_to_stefen_output(basu_output)
output_dict = {}
for key, value in stefen_output.items():
    output_dict.update(dict(zip([mec.dynamicsymbols(key)], [value])))

# Calculation of f_a
B4_num = B4.subs(input_states_dict).subs(para_dict)
B5_num = B5.subs(input_states_dict).subs(para_dict)

# Calculation of dynamic equations
B7_num = B7.subs(input_states_dict).subs(para_dict)
B8_num = B8.subs(T4_zero).subs(input_states_dict).subs(para_dict)

B47 = B4_num.col_join(B7_num)
B58 = B5_num.col_join(B8_num)
output_cal = -B47.inv() * B58
"""


# Steady turning
t5 = time()
# Configuration: lean, steer, pitch
# pitch angle calculated here is to import the function of configuration in the
# steadyturning.py, since the equation it uses to calculate the pitch angle is
# the same as the one from f_c configuration constraint.
lean = pi/8; steer = pi/4
q_conf, q_conf_de = configuration(lean, steer, mp)

# Equilibrium values: u5, T4, u1, u6
u_dep_num_st = u_dep.subs(para_dict).subs(q_conf).subs(steady_conditions)
u_dep_dict_num_st = {udei : u_dep_num_sti[0] for udei, u_dep_num_sti 
                  in zip(u[3:], u_dep_num_st.tolist())}

B8_num_st = B8.subs(steady_conditions).subs(q_conf).subs(para_dict)\
          .subs(qd_kd).subs(u_dep_dict_num_st)

u5sqr_value = sym.solve(B8_num_st[0], u5**2)
if u5sqr_value == []:
    raise ValueError("\nOops! Rear wheel rate in configuration {0} cannot be \
solved. Please select valid configuration according to the plot from <General \
steady turning of a benchmark bicycle model> by Luke. Good luck!\n"\
.format(q_conf))

elif u5sqr_value[0] < 0:
    raise ValueError("Oops! The steady turning in your configuration {0} seems \
Infeasible since no real value appears in rear wheel rate. Please check the \
B8_num_st and try another valid configuraiton according to the plot from \
<General steady turning of a benchmark bicycle model> by Luke.\n"\
.format(q_conf))
else:
    u5_value = -sqrt(u5sqr_value[0])
    print ("\nIt passed feasible steady turning checking and already solved \
the equilibrium values, but it still needs to be checked by eigenvalues of the \ configuration...\n")

equili_dict = {u5: u5_value}

T4_value = sym.solve(B8_num_st[1], T4)[0].subs(equili_dict)

u1_value = u_dep_num_st[0].subs(equili_dict)
u6_value = u_dep_num_st[2].subs(equili_dict)

equili_dict.update(dict(zip([u1, u6, T4], 
                              [u1_value, u6_value, T4_value])))
t6 = time()
print ("Steady turning configuration and equilibrium values generation, \
taking {0} s.".format(t6-t5))

# Contact forces: Fx_r, Fy_r, Fz_r, Fx_f, Fy_f, Fz_f
forces_exp_st = forces.subs(ud_zero).subs(steady_conditions).subs(qd_kd)\
                   .subs(q_conf).subs(equili_dict).subs(para_dict)
force_num_st = sym.solve(forces_exp_st, auxforces)
t7 = time()
print ("Calculating the contact forces for two wheels, taking {0} s. \
".format(t7-t6))

# writing_dyequ(forces)

"""
t8 = time()
# Mass and forcing
# qd_kd, f_a, dynamic  differential about [qd_i, qd_de, ud_i, ud_de]
# qd_kd: M1[qd, ud] = M2
# f_a: M3[qd, ud] = M4
# dynamic: M5[qd, ud] = M6
qud = qd+ud
M1 = sym.zeros(4, len(qud))
M2 = sym.zeros(4, 1)
for i, equi in enumerate(kindiffs):
    for j, qudi in enumerate(qud):
        M1[i, j] = equi.diff(qudi)
    M2[i] = -equi.subs(ud_zero).subs(qd_zero)

M3 = sym.zeros(3, len(qud))
M4 = sym.zeros(3, 1)
for i, equi in enumerate(f_a):
    for j, qudi in enumerate(qud):
        M3[i, j] = equi.diff(qudi)
    M4[i] = -equi.subs(ud_zero).subs(qd_zero)
assert B5 == -M4

M5 = sym.zeros(3, len(qud))
M6 = sym.zeros(3, 1)
for i, equi in enumerate(dynamic):
    for j, qudi in enumerate(qud):
        M5[i, j] = equi.diff(qudi)
    M6[i] = -equi.subs(ud_zero).subs(qd_zero)
assert B8 == -M6

# together
M_left = M1.col_join(M3).col_join(M5)
M_right = M2.col_join(M4).col_join(M6)
# left = M_left.subs(input_states_dict).subs(para_dict)
# right = M_right.subs(input_states_dict).subs(para_dict)
# left.inv()*right #(vs output_dict)
t9 = time()
print ("Generating mass and forcing matrix, taking {0} s".format(t9-t8))


# Simulation
# Import
from scipy.integrate import odeint
from sympy import (Dummy, lambdify)
from bicycle import pitch_from_roll_and_steer
from numpy import linalg
# symbols
dummy_symbols = [Dummy() for i in q+u]
symbols_dict = dict(zip(q+u, dummy_symbols))
symbols_dict.update({T4: 0.})
v = sym.symbols('v')

# Substitution
MM = M_left.subs(symbols_dict).subs(para_dict)
For = M_right.subs(symbols_dict).subs(para_dict)

# Built quick functions for mm and forcing
# Built a function of derivs for q_dot and u_dot
# * here is to convert variables of array (maybe or list) to a tuple for input
# into the lambdify
mm = lambdify(dummy_symbols, MM)
fo = lambdify(dummy_symbols, For)
def derivs(y, t):
    return np.array(linalg.solve(mm(*y), fo(*y))).T[0]

# Initial conditions and time vector defining
# initial condition:
#                  [q1,   q2,   q4,   q3,   u2,   u4,   u5,   u1,   u3,   u6]
#   reference      [0.,   0.,   0.,   lam,  0.,   0.,   v/rr, 0.,   0.,   v/rf]
q1_ref = 0.; q2_ref = 0.; q4_ref = 0.
q3_ref = pitch_from_roll_and_steer(q2_ref, q4_ref, mp['rf'], mp['rr'], 
                                   mp['d1'], mp['d2'], mp['d3'])
q_ini_ref = [q1_ref, q2_ref, q4_ref, q3_ref]

v = 5.0 # m/s
u2_ref = 0.; u4_ref = 0.; u5_ref = v/mp['rr']
u1_ref = 0.; u3_ref = 0.; u6_ref = v/mp['rf']
u_ini_ref = [u2_ref, u4_ref, u5_ref, u1_ref, u3_ref, u6_ref]

y0_ref = np.array(q_ini_ref + u_ini_ref)
t = np.linspace(0, 5, 1000)

# Intgration: y: (10, 1)
# reorder the y to correspond [q1, q2, q3, q4, u1, u2, u3, u4, u5, u6]
#                  instead of [q1, q2, q4, q3, u2, u4, u5, u1, u3, u6]
y_ref = odeint(derivs, y0_ref, t)
y_ref_reorder = np.hstack((y_ref[:, :2], y_ref[:,3:4], y_ref[:,2:3], 
                  y_ref[:,5:6], y_ref[:,7:9], y_ref[:,4:5], y_ref[:,6:7], 
                  y_ref[:,9:]))

# Plotting
from matplotlib.pyplot import (figure, plot, legend, xlabel, ylabel, show,
                                title)
for i in range(len(q)):
    figure(1)
    plot(t, y_ref_reorder[:, i], label='q'+str(i+1))
for j in range(len(u)):
    figure(2)
    plot(t, y_ref_reorder[:, j+len(q)], label='u'+str(j+1))

figure(1)
title('Angle performance in reference configuration')
legend(loc=0)
xlabel('Time (s)')
ylabel('Angle (rad)')

figure(2)
title('Rate performance in reference configuration')
legend(loc=0)
xlabel('Time (s)')
ylabel('Rate (rad/s)')
show()

t10 = time()
print ('Simulation and plot, taking {0} s'.format(t10-t9))
"""
