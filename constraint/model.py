import sympy as sym
from numpy import *
import sympy.physics.mechanics as mec
import pdb

import os

from sympy import signsimp, factor_terms

mec.Vector.simp = False
mec.mechanics_printing()


class bicycle_model:
    def __init__(self):

        #=================
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


        #===================================
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

        self.coordinatesInde = [q1, q2, q4]
        self.coordinatesDe = [q3]
        self.speedsInde = [u2, u4, u5]
        self.speedsDe = [u1, u3, u6]

        #==================================
        # Axiliary speeds at contact points

        #rear wheel
        ua1, ua2 = mec.dynamicsymbols ('ua1 ua2')
        #front wheel
        ua4, ua5 = mec.dynamicsymbols ('ua4 ua5')

        self.auxiliarySpeeds = [ua1, ua2, ua4, ua5]

        #================================
        # Orientation of Reference Frames

        # bicycle frame yaw
        A.orient(N, 'Axis', [q1, N['3']])
        # bicycle frame roll
        B.orient(A, 'Axis', [q2, A['1']])
        # bicycle frame pitch
        C.orient(B, 'Axis', [q3, B['2']])
        # fork/handlebar steer
        E.orient(C, 'Axis', [q4, C['3']])


        #===========
        # parameters

        # geometry
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

        rf, rr = sym.symbols('rf rr')
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


        #===================
        # Angular Velocities

        A.set_ang_vel(N, u1 * N['3'])
        B.set_ang_vel(A, u2 * A['1'])
        C.set_ang_vel(B, u3 * B['2'])
        E.set_ang_vel(C, u4 * C['3'])

        D.set_ang_vel(C, u5 * C['2'])
        F.set_ang_vel(E, u6 * E['2'])

        #=====================
        # special unit vectors

        #pitch back for the unit vector g_3 along front wheel radius;
        g_3 =  (mec.express(A['3'], E) - mec.dot(E['2'], A['3'])*E['2']).normalize() 
        #another way: g_3 = E['2'].cross(A['3']).cross(E['2']).normalize() 


        #roll back for longitudinal and lateral unit vector of front wheel
        long_v = mec.cross (E['2'], A['3']).normalize()
        lateral_v = mec.cross (A['3'], long_v).normalize() 


        #======================
        # points and velocities

        #rear wheel contact point, dn
        dn = mec.Point('dn')
        dn.set_vel(N, ua1 * A['1'] + ua2 * A['2']) 
        #dn.set_vel(N, ua1 * N['1'] + ua2 * N['2'])


        # rear wheel center, do
        do = dn.locatenew('do', -rr * B['3'])
        #do = dn.locatenew('do', -rtr * A['3'] - rr * B['3']) 
        #rear wheel center with rear tire radius rtr

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


        #Front wheel contact point, fn
        fn = mec.Point('fn')
        fn.set_vel(N, ua4 * long_v + ua5 * lateral_v)
        #fn.set_vel(N, ua4 * N['1'] + ua5 * N['2'])


        #front wheel center
        fo = fn.locatenew('fo', -rf * g_3)
        #fo = fn.locatenew('fo', -ftr * A['3'] - rf * g_3) 
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

        #=============================
        # Holo and nonholo Constraints

        self.holonomic = [fn.pos_from(dn).dot(A['3'])]

        self.nonholonomic = [(v_SAF_1-v_SAF_2).dot(uv) for uv in E]


        #=============
        # Rigid Bodies

        #Inertia
        Ic = mec.inertia(C, ic11, ic22, ic33, 0.0, 0.0, ic31)
        Id = mec.inertia(C, id11, id22, id11, 0.0, 0.0, 0.0) #rear wheel
        Ie = mec.inertia(E, ie11, ie22, ie33, 0.0, 0.0, ie31)
        If = mec.inertia(E, if11, if22, if11, 0.0, 0.0, 0.0) #front wheel

        #bodies
        rearFrame_inertia = (Ic, co)
        rearFrame=mec.RigidBody('rearFrame',co,C,mc,rearFrame_inertia)

        rearWheel_inertia = (Id, do)
        rearWheel=mec.RigidBody('rearWheel',do,D,md,rearWheel_inertia)

        frontFrame_inertia = (Ie, eo)
        frontFrame=mec.RigidBody('frontFrame',eo,E,me,frontFrame_inertia)

        frontWheel_inertia = (If, fo)
        frontWheel=mec.RigidBody('frontWheel',fo,F,mf,frontWheel_inertia)

        bodyList = [rearFrame, rearWheel, frontFrame, frontWheel]



        #==========================
        # Generalized Active Forces

        # steer torque
        T4 = mec.dynamicsymbols('T4')

        #road_wheel contact forces
        Fx_r, Fy_r, Fx_f, Fy_f = mec.dynamicsymbols('Fx_r Fy_r Fx_f Fy_f')


        # gravity
        Fco = (co, mc * g * A['3'])
        Fdo = (do, md * g * A['3'])
        Feo = (eo, me * g * A['3'])
        Ffo = (fo, mf * g * A['3'])

        #road_wheel contact forces
        F_r = (dn, Fx_r * A['1'] + Fy_r * A['2'])
        F_f = (fn, Fx_f * long_v + Fy_f * lateral_v)

        #F_r = (dn, Fx_r * N['1'] + Fy_r * N['2'])
        #F_f = (fn, Fx_f * N['1'] + Fy_f * N['2'])

        # input torques
        Tc = (C, - T4 * C['3']) #back steer torque to rear frame

        Te = (E, T4 * C['3']) #steer torque to front frame


        forceList = [Fco, Fdo, Feo, Ffo, F_r, F_f, Tc, Te]


        #===================================    
        # Kinematical Differential Equations

        kinematical = [q1d - u1,
                       q2d - u2,
                       q3d - u3,
                       q4d - u4]

        #=======================
        print ('Kanes Method')

        self.kane = mec.KanesMethod(N, q_ind=self.coordinatesInde, u_ind=self.speedsInde, 
            kd_eqs=kinematical, q_dependent=self.coordinatesDe, 
            configuration_constraints=self.holonomic, 
            u_dependent=self.speedsDe, velocity_constraints=self.nonholonomic, 
            u_auxiliary=self.auxiliarySpeeds)
            
        #reminder: u1--yaw rate, u2--roll rate, u3--pitch rate, u4--steer rate, 
        #u5--rear wheel ang. vel., u6--front wheel ang. vel.

        (fr, frstar)= self.kane.kanes_equations(forceList, bodyList)

        self.kdd = self.kane.kindiffdict()


    def mass_matrix_full(self):
        self.mmFull = self.kane.mass_matrix_full.subs(self.kdd)


    def forcing_full(self):
        self.forceFull = self.kane.forcing_full.subs(self.kdd)


    def linearized_a(self):
        """linearization of focing matrix, obtaining A matrix."""
        self.forcingLinA = self.kane.linearize()[0].subs(self.kdd)


    def linearized_b(self):
        """linearization of forcing matrix, obtaining  B matrix."""
        self.forcingLinB = self.kane.linearize()[1].subs(self.kdd)

    def parameters_symbols(self, mooreParameters):
        """Returns a dictionary of parameters whose keys are symbols instead of
        strings in mooreParameters."""
        mp = mooreParameters
        self.parameters = {}
        for key, value in mp.items():
            self.parameters.update(dict(zip([sym.symbols(key)], [value])))

    def coordinates_dynamicsymbols(self, coordinates):
        """Returns a dictionary of coordinates whose keys are dynamic symbols 
        instead of strings of coordinates."""
        self.coordinates = {}
        for key, value in coordinates.items():
            self.coordinates.update(dict(zip([mec.dynamicsymbols(key)], [value])))

    def speeds_dynamicsymbols(self, speeds):
        """Returns a dictionary of speeds whose keys are dynamic symbols 
        instead of strings of speeds."""
        self.speeds = {}
        for key, value in speeds.items():
            self.speeds.update(dict(zip([mec.dynamicsymbols(key)], [value])))

    def auxiliary_speeds_zero(self):
        """Returns a dictionary of zero auxiliary speeds."""
        self.auxiliarySpeedsZeros = dict(zip(self.auxiliarySpeeds, 
                                            zeros(len(self.auxiliarySpeeds))))

    def linearized_reference_configuration(self, lam, rR, rF):
        """Returns the linearized model at the reference configuration."""
        v = sym.Symbol('v')

        q1, q2, q4 = self.coordinatesInde
        q3 = self.coordinatesDe[0]
        u2, u4, u5 = self.speedsInde
        u1, u3, u6 = self.speedsDe

        self.referenceConfiguration = {q1: 0., q2: 0., q3:lam, q4:0., u1:0., 
        u2:0., u3:0., u4:0., u5: -v/rR, u6: -v/rF}

