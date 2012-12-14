#========================
#print ('built functions')

def configuration(lean, steer):
    """Returns the configuration of steady turning in lean, steer
    and pitch angle in radians and degrees by giving front twos values.

    Return
    ------
    q_dict: dictionary
        The dictionary of steady turning configuration.
    q_dict_d: dictionary
        The degrees units of q_dict.
    """
    pitch = bi.pitch_from_roll_and_steer(lean, steer, mp['rf'], mp['rr'], 
            mp['d1'], mp['d2'], mp['d3'])
    
    q_dict = {q2: lean, q3:pitch, q4:steer}
    
    lean_d, pitch_d, steer_d = [value*180/pi for value in [lean, pitch, steer]]
    
    q_dict_d = {q2: lean_d, q3: pitch_d, q4: steer_d}
    
    return q_dict, q_dict_d

def forcing_dynamic_equations(forcingEquations, parameters, qDict, uDict):
    """Returns the focing dynamic equations in specific configuration.
    Note
    ----
    There are no u_dots here.

    Parameter
    ---------
    forcingEquatioins: dictionary
        Forcing matrix eqautions.
    Parameters: dictionary
        Bicycle parameters.
    qDict: dictionary
        Configuration of steady turning: lean, steer and pitch angles.
    uDict: dictionary
        zeros u's: lean rate, pitch rate, steer rate.

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
    uDict: dictionary
        A dictionary consists of zeros u values for steady turning, 
        leanrate(u2), pitchrate(u3), and steerrate(u4).

    Return
    ------
    inde_expression: array, 3-by-1
        A list of expressions associated with all independent generalized speeds.
        The expressions correspond to all dependent speeds.
    inde_expression_subs: array, 1-by-3
        Substituding the uDict into the inde_expression.
    """
    #independent generalized speeds
    nonho_coeff_inde_sym = [[value.expand().coeff(u2), value.expand().coeff(u4), 
                             value.expand().coeff(u5)] for value in nonholonomic]
    #dependent generalized speeds
    nonho_coeff_de_sym = [[value.expand().coeff(u1), value.expand().coeff(u3), 
                           value.expand().coeff(u6)] for value in nonholonomic]
    
    nonho_coeff_inde_value = [[value.subs(parameters).subs(qDict) for value in value_1] 
                              for value_1 in nonho_coeff_inde_sym]
    
    nonho_coeff_de_value = [[value.subs(parameters).subs(qDict) for value in value_1] 
                            for value_1 in nonho_coeff_de_sym]
    
    #independent speeds
    inde_states = matrix([[u2], [u4], [u5]])
    
    nonho_coeff_inde_ma = asmatrix(nonho_coeff_inde_value)
    nonho_coeff_de_ma = asmatrix(nonho_coeff_de_value)

    inde_expression_ma = (nonho_coeff_de_ma.I) * (- nonho_coeff_inde_ma * inde_states)
    
    inde_expression = asarray(inde_expression_ma)
    
    inde_expression_list = [value[0].subs(u_dict) for value in inde_expression]
    
    return inde_expression, inde_expression_list

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
