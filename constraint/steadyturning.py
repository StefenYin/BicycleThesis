import bicycle as bi

import sympy as sym
import sympy.physics.mechanics as mec
from numpy import *


q2, q3, q4 = mec.dynamicsymbols('q2 q3 q4')
u1, u2, u3, u4, u5, u6= mec.dynamicsymbols('u1 u2 u3 u4 u5 u6')
T4 = mec.dynamicsymbols('T4')


def configuration(lean, steer, mooreParameters):
    """Returns the configuration of steady turning in lean, steer and pitch 
    angle in radians and degrees by giving front twos values and bicycle
    parameters in Jason's set.

    Parameter
    ---------
    lean: float
        A given lean angle.
    steer: float
        A given steer angle.
    mooreParameters: dictionary
        Bicycle parameters in Moore's set, but keys are string.

    Return
    ------
    q_dict: dictionary
        The dictionary of steady turning configuration.
    q_dict_d: dictionary
        The degrees units of q_dict.

    """

    mp = mooreParameters

    pitch = bi.pitch_from_roll_and_steer(lean, steer, mp['rf'], mp['rr'], 
            mp['d1'], mp['d2'], mp['d3'])

    lean_d, pitch_d, steer_d = [value*180/pi for value in [lean, pitch, steer]]

    q_dict = {q2: lean, q3:pitch, q4:steer}
    q_dict_d = {q2: lean_d, q3: pitch_d, q4: steer_d}
    
    return q_dict, q_dict_d


def speeds_zeros(dynamicSymbolsSpeeds):
    """Returns the zeros of speeds, e.g. lean rate, pitch rate, and steer rate 
    in steady turning.

    Parameter
    ---------
    dynamicSymbolsSpeeds: list
        A list of dynamic symbols of speeds.

    Return
    ------
    u_zeros_dict : dictionary
        Zero speeds.

    """

    u_zeros_dict = dict(zip(dynamicSymbolsSpeeds,zeros(len(dynamicSymbolsSpeeds))))

    return u_zeros_dict


def forcing_dynamic_equations(forcingEquations, parameters, qDict, uDict):
    """Returns the focing dynamic equations in specific configuration.

    Note
    ----
    There are no u_dots here.

    Parameter
    ---------
    forcingEquatioins: list
        Forcing matrix equations derived from Kane' Method with Sympy.
    parameters: dictionary
        Bicycle parameters in Moore' set, but keys are symbols in forcingEquations.
    qDict: dictionary
        Configuration of steady turning: lean, steer and pitch angles.
    uDict: dictionary
        zeros u's of steady turning: lean rate, pitch rate, steer rate.

    Return
    ------
    dynamicEquations: dictionary
        Dynamic equations with yaw rate, rear wheel rate, front wheel rate, and 
        steer torque.

    """

    F_full = forcingEquations.subs(uDict).subs(qDict).subs(parameters).expand()
    
    dynamicEquations = F_full[4:7]
    
    return dynamicEquations


def dynamic_equations_coefficients(dynamicEquations):
    """Returns the coeffcients in the following form:
    a*(u1*u1) + b*(u1*u5) + c*(u1*u6) + d*T4 + e = 0.
    
    Parameter
    ---------
    dynamicEquations: dictionary
        Dynamic equations in above form after expand.

    Return
    ------
    coefficients: dictionary
        a, b, c, d, e in a row for each dynamic equation.

    """

    num_dict = {u1: 0.0, u5: 0.0, u6: 0.0, T4: 0.0}
    
    dynamicEquations = [value.expand() for value in dynamicEquations]
    
    coefficients = [[value.coeff(u1*u1), value.coeff(u1*u5), 
                     value.coeff(u1*u6), value.coeff(T4), 
                     value.subs(num_dict)] for value in dynamic_equ]
    
    return coefficients


def de_by_inde(nonholonomic, qDict, parameters, uDict):
    """Returns the dependent generalized speeds (u1, u3, u6) expressed by 
    independent speeds (u2, u4, u5).
    
    Parameter
    ---------
    nonholonomic: array-like
        A list of nonholonomic constraint equations.

    Return
    ------
    inde_expression: array, 3-by-1
        A list of expressions associated with all independent generalized speeds.
        The expressions correspond to all dependent speeds.
    inde_expression_list: a list
        Substituding the uDict into the inde_expression.

    """

    # Independent generalized speeds coefficients: u2, u4 and u5
    # Dependent generalized speeds coefficients: u1, u3 and u6
    nonho_coeff_inde_sym = [[value.expand().coeff(u2), value.expand().coeff(u4), 
                             value.expand().coeff(u5)] for value in nonholonomic]

    nonho_coeff_de_sym = [[value.expand().coeff(u1), value.expand().coeff(u3), 
                           value.expand().coeff(u6)] for value in nonholonomic]

    nonho_coeff_inde_value = [[value.subs(parameters).subs(qDict) for value in value_1] 
                              for value_1 in nonho_coeff_inde_sym]

    nonho_coeff_de_value = [[value.subs(parameters).subs(qDict) for value in value_1] 
                            for value_1 in nonho_coeff_de_sym]

    # Constraint equations -> de expressed by inde
    inde_states = matrix([[u2], [u4], [u5]])

    nonho_coeff_inde_ma = asmatrix(nonho_coeff_inde_value)
    nonho_coeff_de_ma = asmatrix(nonho_coeff_de_value)

    inde_expression_ma = (nonho_coeff_de_ma.I) * (- nonho_coeff_inde_ma * inde_states)
    inde_expression_arr = asarray(inde_expression_ma)

    inde_expression = [value[0] for value in inde_expression_arr]
    inde_expression_subs = [value[0].subs(uDict) for value in inde_expression_arr]

    return inde_expression, inde_expression_subs


def dynamic_nonholonomic_equations(indeExpression_list, dynamicEquations):
    """Return the final equations after substitude independent expressions 
    into the dynamic equations.

    Note
    ----
    The equations of indeExpression and dynamicEquations are already plugged 
    with parameters, qDict and uDict.

    Return
    ------
    dynamic_nonho_equ: array-like
        Dynamic equations after substituding uDict and indeExpression.

    """

    inde_dict = dict(zip([u1, u3, u6], indeExpression_list))
    
    dynamic_nonho_equ = [value.subs(inde_dict).expand() for value in dynamicEquations]

    return dynamic_nonho_equ
