"""
Model module:
1, building model;
2, mass_matrix_full, forcing_full, linearized;
3, parameters, coordinates, speeds convert from string to symbols;
4, zero auxiliary speeds, reference configuration;
5, contact forces.
"""

import os
import pdb

import sympy as sym
import sympy.physics.mechanics as mec
from numpy import *
from sympy import (signsimp, factor_terms)

mec.Vector.simp = False
mec.mechanics_printing()



class BicycleModel(object):
    """Create a whipple model with auxiliary speeds on the contact points."""


    def __init__(self):

        # Generalized Coordinates and Speeds:
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

        self._coordinatesInde = [q1, q2, q4]
        self._coordinatesDe = [q3]
        self._coordinates = [q1, q2, q3, q4]
        self._speedsInde = [u2, u4, u5]
        self._speedsDe = [u1, u3, u6]
        self._speeds = [u1, u2, u3, u4, u5, u6]
        self._speedsDerivative = [u1d, u2d, u3d, u4d, u5d, u6d]

        # Axiliary speeds at contact points:
        # Rear wheel: ua1, ua2
        # Front wheel: ua4, ua5
        ua1, ua2, ua3 = mec.dynamicsymbols ('ua1 ua2 ua3')
        ua4, ua5, ua6 = mec.dynamicsymbols ('ua4 ua5 ua6')
        ua1d, ua2d, ua3d = mec.dynamicsymbols ('ua1 ua2 ua3',1)
        ua4d, ua5d, ua6d = mec.dynamicsymbols ('ua4 ua5 ua6',1)

        self._auxiliarySpeeds = [ua1, ua2, ua3, ua4, ua5, ua6]
        self._auxiliarySpeedsDerivative = [ua1d, ua2d, ua3d, ua4d, ua5d, ua6d]
        ua_zero = {uai: 0 for uai in self._auxiliarySpeeds}

        # Reference Frames:
        # Newtonian Frame: N
        # Yaw Frame: A
        # Roll Frame: B
        # Pitch & Bicycle Frame: C
        # Steer & Fork/Handlebar Frame: E
        # Rear Wheel Frame: D
        # Front Wheel Frame: F
        N = mec.ReferenceFrame('N', indices=('1', '2', '3'))
        A = mec.ReferenceFrame('A', indices=('1', '2', '3'))
        B = mec.ReferenceFrame('B', indices=('1', '2', '3'))
        C = mec.ReferenceFrame('C', indices=('1', '2', '3'))
        E = mec.ReferenceFrame('E', indices=('1', '2', '3'))
        D = mec.ReferenceFrame('D', indices=('1', '2', '3'))
        F = mec.ReferenceFrame('F', indices=('1', '2', '3'))

        # Orientation of Reference Frames:
        # bicycle frame yaw: N->A
        # bicycle frame roll: A->B
        # pitch to rear frame: B->C
        # fork/handlebar steer: C->E
        A.orient(N, 'Axis', [q1, N['3']])
        B.orient(A, 'Axis', [q2, A['1']])
        C.orient(B, 'Axis', [q3, B['2']])
        E.orient(C, 'Axis', [q4, C['3']])

        # Angular Velocities and define the generalized speeds:
        A.set_ang_vel(N, u1 * N['3'])
        B.set_ang_vel(A, u2 * A['1'])
        C.set_ang_vel(B, u3 * B['2'])
        E.set_ang_vel(C, u4 * C['3'])
        D.set_ang_vel(C, u5 * C['2'])
        F.set_ang_vel(E, u6 * E['2'])

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
        self._massSym = [mc, md, me, mf]

        ic11, ic22, ic33, ic31 = sym.symbols('ic11 ic22 ic33 ic31') #rear frame
        id11, id22 = sym.symbols('id11 id22')  #rear wheel
        ie11, ie22, ie33, ie31 = sym.symbols('ie11 ie22 ie33 ie31')  #front frame
        if11, if22 = sym.symbols('if11 if22') #front wheel

        # Special unit vectors:
        # g_3: direction along front wheel radius.
        # long_v: longitudinal direction of front wheel.
        # lateral_v: lateral direction of front wheel.
        g_3 =  (mec.express(A['3'], E) - mec.dot(E['2'], A['3'])*E['2']).normalize() 
        # Or g_3 = E['2'].cross(A['3']).cross(E['2']).normalize() 
        long_v = mec.cross(E['2'], g_3).normalize()
        # E.y ^ g_3 # ^ -> cross,  & -> dot
        lateral_v = mec.cross(A['3'], long_v)

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

        # Rear
        dn = mec.Point('dn')
        dn.set_vel(N, ua1 * A['1'] + ua2 * A['2'] + ua3 * A['3']) 

        do = dn.locatenew('do', -rr * B['3'])
        do.v2pt_theory(dn, N, D)
        do.set_acc(N, do.vel(N).subs(ua_zero).diff(t, C) + 
           mec.cross(C.ang_vel_in(N), do.vel(N).subs(ua_zero)))

        co = mec.Point('co')
        co.set_pos(do, l1 * C['1'] + l2 * C['3'])
        co.v2pt_theory(do, N, C)
        co.set_acc(N, co.vel(N).subs(ua_zero).diff(t, C) + 
           mec.cross(C.ang_vel_in(N), co.vel(N).subs(ua_zero)))

        # Front
        fn = mec.Point('fn')
        fn.set_vel(N, ua4 * long_v + ua5 * lateral_v + ua6 * A['3'])

        fo = fn.locatenew('fo', -rf * g_3)
        # OR fo = fn.locatenew('fo', -ftr * A['3'] - rf * g_3) 
        fo.v2pt_theory(fn, N, F)
        fo.set_acc(N, fo.vel(N).subs(ua_zero).diff(t, E) + 
           mec.cross(E.ang_vel_in(N), fo.vel(N).subs(ua_zero)))

        eo = mec.Point('eo')
        eo.set_pos(fo, l3 * E['1'] + l4 * E['3'])
        eo.v2pt_theory(fo, N, E)
        eo.set_acc(N, eo.vel(N).subs(ua_zero).diff(t, E) + 
           mec.cross(E.ang_vel_in(N), eo.vel(N).subs(ua_zero)))

        SAF = do.locatenew('SAF', d1 * C['1'] + d2 * E['3'])
        SAF.set_pos(fo, -d3 * E['1'])

        # Velociy of SAF in two ways:
        # OR v_SAF_1 = v2pt_theory(do, N, C)
        # OR v_SAF_2 = v2pt_theory(fo, N, E)
        v_SAF_1 = do.vel(N) + mec.cross(C.ang_vel_in(N), SAF.pos_from(do))
        v_SAF_2 = fo.vel(N) + mec.cross(E.ang_vel_in(N), SAF.pos_from(fo))

        # Holo and nonholo Constraints:
        self._holonomic = [fn.pos_from(dn).dot(A['3'])]
        self._nonholonomic = [(v_SAF_1-v_SAF_2).dot(uv) for uv in E]

        # Rigid Bodies:
        # Inertia: Ic, Id, Ie, If
        # Bodies: rearFrame, rearWheel, frontFrame, frontWheel
        Ic = mec.inertia(C, ic11, ic22, ic33, 0.0, 0.0, ic31)
        Id = mec.inertia(C, id11, id22, id11, 0.0, 0.0, 0.0)
        Ie = mec.inertia(E, ie11, ie22, ie33, 0.0, 0.0, ie31)
        If = mec.inertia(E, if11, if22, if11, 0.0, 0.0, 0.0)

        rearFrame_inertia = (Ic, co)
        rearFrame=mec.RigidBody('rearFrame',co,C,mc,rearFrame_inertia)
        rearWheel_inertia = (Id, do)
        rearWheel=mec.RigidBody('rearWheel',do,D,md,rearWheel_inertia)
        frontFrame_inertia = (Ie, eo)
        frontFrame=mec.RigidBody('frontFrame',eo,E,me,frontFrame_inertia)
        frontWheel_inertia = (If, fo)
        frontWheel=mec.RigidBody('frontWheel',fo,F,mf,frontWheel_inertia)

        bodyList = [rearFrame, rearWheel, frontFrame, frontWheel]

        # Generalized Active Forces:
        # T4: steer torque.
        # Fx_r, Fy_r: rear wheel contact forces.
        # Fx_f, Fy_f: front wheel contact forces.
        # Fco, Fdo, Feo, Ffo: gravity of four rigid bodies of bicycle.
        T4 = mec.dynamicsymbols('T4')
        Fx_r, Fy_r, Fz_r, Fx_f, Fy_f, Fz_f= mec.dynamicsymbols('Fx_r Fy_r Fz_r Fx_f Fy_f Fz_f')

        Tc = (C, - T4 * C['3']) #back steer torque to rear frame
        Te = (E, T4 * C['3']) #steer torque to front frame

        F_r = (dn, Fx_r * A['1'] + Fy_r * A['2'] + Fz_r * A['3'])
        F_f = (fn, Fx_f * long_v + Fy_f * lateral_v + Fz_f * A['3'])
        # OR F_r = (dn, Fx_r * N['1'] + Fy_r * N['2'])
        # OR F_f = (fn, Fx_f * N['1'] + Fy_f * N['2'])

        Fco = (co, mc * g * A['3'])
        Fdo = (do, md * g * A['3'])
        Feo = (eo, me * g * A['3'])
        Ffo = (fo, mf * g * A['3'])

        forceList = [Fco, Fdo, Feo, Ffo, F_r, F_f, Tc, Te]

        self._inputForces = [T4]
        self._auxiliaryForces = [Fx_r, Fy_r, Fz_r, Fx_f, Fy_f, Fz_f]

        # Kinematical Differential Equations:
        kinematical = [q1d - u1,
                       q2d - u2,
                       q3d - u3,
                       q4d - u4]

        # Kanes Method:
        self._kane = mec.KanesMethod(
            N, q_ind=self._coordinatesInde, u_ind=self._speedsInde, kd_eqs=kinematical, 
            q_dependent=self._coordinatesDe, configuration_constraints=self._holonomic, 
            u_dependent=self._speedsDe, velocity_constraints=self._nonholonomic, 
            u_auxiliary=self._auxiliarySpeeds
            )

        (fr, frstar)= self._kane.kanes_equations(forceList, bodyList)

        self._Fr = fr
        self._Fr_star = frstar
        self._kdd = self._kane.kindiffdict()

        # Setting another contact points for calculating turning radius:
        # Turning radius: rear wheel: Rr, front wheel Rf;
        # Contact points: rear contact: P; front contact: Q; Turning center: TC;
        # Relative position: P_Q_A1, P_Q_A2.
        Rr, Rf = sym.symbols('Rr Rf')
        P = mec.Point('P')
        TC = P.locatenew('TC', Rr*A['2'])
        Q = TC.locatenew('Q', -Rf*lateral_v)
        self._turningRadiusSym = [Rr, Rf]

        P_Q_A1 = Q.pos_from(P).dot(A['1'])
        P_Q_A2 = Q.pos_from(P).dot(A['2'])
        self._contact_posi_pq = [P_Q_A1, P_Q_A2]

        fn_dn_A1 = fn.pos_from(dn).dot(A['1'])
        fn_dn_A2 = fn.pos_from(dn).dot(A['2'])
        self._contact_posi_dfn = [fn_dn_A1, fn_dn_A2]

        # Total center of mass:
        # individual center of mass of four rigid bodies
        co_dn_A = [mec.dot(co.pos_from(dn), A_x) for A_x in A]
        do_dn_A = [mec.dot(do.pos_from(dn), A_x) for A_x in A]
        fo_dn_A = [mec.dot(fo.pos_from(dn), A_x) for A_x in A]
        eo_dn_A = [mec.dot(eo.pos_from(dn), A_x) for A_x in A]

        self._bodies_dn_A = array(
                            [[co_dn_A[0], do_dn_A[0], eo_dn_A[0], fo_dn_A[0]], 
                            [co_dn_A[1], do_dn_A[1], eo_dn_A[1], fo_dn_A[1]],
                            [co_dn_A[2], do_dn_A[2], eo_dn_A[2], fo_dn_A[2]]]
                            )


    def mass_matrix_full(self):
        """Returns mass matrix."""

        self.mmFull = self._kane.mass_matrix_full.subs(self._kdd)

        return self.mmFull


    def forcing_full(self):
        """Returns forcing matrix."""

        self.forceFull = self._kane.forcing_full.subs(self._kdd)

        return self.forceFull


    def linearized_a(self):
        """linearization of focing matrix, obtaining A matrix."""

        self.forcingLinA = self._kane.linearize()[0].subs(self._kdd)

        return self.forcingLinA


    def linearized_b(self):
        """linearization of forcing matrix, obtaining  B matrix."""

        self.forcingLinB = self._kane.linearize()[1].subs(self._kdd)

        return self.forcingLinB


    def linearized_reference_configuration(self, lam, rR, rF):
        """Returns the linearized model at the reference configuration."""

        v = sym.Symbol('v')

        q1, q2, q4 = self._coordinatesInde
        q3 = self._coordinatesDe[0]
        u2, u4, u5 = self._speedsInde
        u1, u3, u6 = self._speedsDe

        self.referenceConfiguration = {q1: 0., q2: 0., q3:lam, q4:0., u1:0., 
        u2:0., u3:0., u4:0., u5: -v/rR, u6: -v/rF}

        return self.referenceConfiguration


    def contact_forces(self):
        """Returns contact forces on each wheel."""

        self.conForceNoncontri = self._kane.auxiliary_eqs.applyfunc(
                            lambda w: factor_terms(signsimp(w))).subs(self._kdd)

        return self.conForceNoncontri



def strings2symbols(strings, go2type = None):
    """Returns a dictionary with keys being symbols instead of strings.
    Symbols here can be ordinary symbols or dynamic symbols.

    Parameter
    ---------
    strings: a dictionary
        A dictionary with keys being strings.
    go2type: string
        A type of keys of the dictionary to be transfered into.
        Two options only: "orsymbols" for ordinary symbols, 
                          "dysymbols" for dynamic symbols.

    Return
    ------
    symbols: a dictionay
        A dictionary with keys being symbols.

    """

    symbols = {}

    if go2type is None:
        raise TypeError("The type of keys in the dictionary to be transfered "
                        "into must supply")

    elif go2type is "orsymbols":
        for key, value in strings.items():
            symbols.update(dict(zip([sym.symbols(key)], [value])))

    elif go2type is "dysymbols":
        for key, value in strings.items():
            symbols.update(dict(zip([mec.dynamicsymbols(key)], [value])))

    return symbols


def zeros_dict(symbols):
    """Returns a dictionary with keys being symbols, values being zeros.

    Parameters
    ----------
    symbols: a list
        A list of symbols or strings.

    """

    zeros_symbols = dict(zip(symbols, zeros(len(symbols))))

    return zeros_symbols
