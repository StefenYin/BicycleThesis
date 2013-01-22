"""
In this model, there is one function for writing the dynamic equations from
kane's method into a file.
"""

import textwrap
from sympy.physics.mechanics import dynamicsymbols
import os

def writing_dyequ(dyEquations, path = None, symbols = None):
    """Return a file containing dynamic equations following the provided path 
    (or default path).

    Parameters
    ----------
    dyEquations : a list
        A list of dynamic equations derived from kane' class in mechanics which
        is ready to be written into the path.
    path : a string
        A string of path. The default path is going to be the path in my
        computers.
    symbols : a list
        A list strings of symbols which is going to be transfered into the
        dynamic symbols.
        Note: the symbols here is extra symbols except common q's, u's defined
        in the model. e.g. maybe "at" for acceleration of total center of mass.
        But, here in the function, most of symbols will be covered.

    """

    equ = dyEquations

    # Path
    if path is None:
        path1 = '/home/stefenyin'
        path2 = '/home/stefenstudy'
        if os.path.exists(path1):
            path = path1
        else:
            path = path2
        path = path + '/bicycle/bi_equations_writing/writing.txt'
    else:
        path = path

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

    # Strings of ordinary symbols and dynamic symbols
    if symbols is None:
        symbols_orstr = []
        symbols_dystr = []
    else:
        symbols_orstr = symbols
        symbols_dystr = [str(dynamicsymbols(symi)) for symi in symbols_orstr]
    
    ud_dystr = ['u2d(t)','u4d(t)','u5d(t)','u1d(t)','u3d(t)','u6d(t)']
    qu_dystr = ['q1(t)', 'q2(t)', 'q4(t)', 'q3(t)',
              'u2(t)', 'u4(t)', 'u5(t)', 'u1(t)', 'u3(t)', 'u6(t)']
    forces_dystr = ['Fx_r(t)', 'Fy_r(t)',  'Fz_r(t)',
                    'Fx_f(t)', 'Fy_f(t)',  'Fz_f(t)',
                    'T4(t)']

    ud_orstr = ['u2d', 'u4d', 'u5d', 'u1d', 'u3d', 'u6d']
    qu_orstr = ['q1', 'q2', 'q4', 'q3',
                  'u2', 'u4', 'u5', 'u1', 'u3', 'u6']
    forces_orstr = ['Fx_r', 'Fy_r',  'Fz_r',
                    'Fx_f', 'Fy_f',  'Fz_f',
                    'T4']

    orstr = symbols_orstr + ud_orstr + qu_orstr + forces_orstr
    dystr = symbols_dystr + ud_dystr + qu_dystr + forces_dystr

    para_symbols = ['rF','rR', 
                    'd1','d2','d3', 
                    'l1','l2','l3', 'l4',
                    'g', 
                    'mc','md','me','mf', 
                    'ic11','ic22','ic33','ic31', 
                    'id11','id22', 
                    'ie11','ie22','ie33','ie31', 
                    'if11','if22',
                    'v']

    # Wrap building
    wrapper = textwrap.TextWrapper(width=73, initial_indent='   ',
                                    subsequent_indent='        ')

    # Start writing
    # Write Head
    f.write('Here are the dynamic equations.\n\n')

    for equi in equ:
        exp_str = str(equi.evalf(n=3))

        title = 'In equation: '
        para = 'Parameters: '
        signals = 'Dynamic symbols: '

        for i, j in zip(orstr, dystr):
            if j in exp_str:
                signals += j + ', '
            exp_str = exp_str.replace(j, i)

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
        f.write (note + '.\n')

        exp_str = wrapper.wrap(exp_str)
        for j in range(len(exp_str)):
            f.write (exp_str[j] + ' \\' + '\n')

        f.write ('\n\n')

    f.close()
