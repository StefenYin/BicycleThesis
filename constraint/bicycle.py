from numpy import pi, sin, cos, tan
import numpy as np

def benchmark_parameters():
    """Returns the benchmark bicycle parameters from [Meijaard2007]_."""

    p = {}

    p['w'] = 1.02
    p['c'] = 0.08
    p['lam'], p['lambda'] = pi / 10., pi / 10.
    p['g'] = 9.81
    p['rR'] = 0.3
    p['mR'] = 2.0
    p['IRxx'] = 0.0603
    p['IRyy'] = 0.12
    p['xB'] = 0.3
    p['zB'] = -0.9
    p['mB'] = 85.0
    p['IBxx'] = 9.2
    p['IByy'] = 11.0
    p['IBzz'] = 2.8
    p['IBxz'] = 2.4
    p['xH'] = 0.9
    p['zH'] = -0.7
    p['mH'] = 4.0
    p['IHxx'] = 0.05892
    p['IHyy'] = 0.06
    p['IHzz'] = 0.00708
    p['IHxz'] = -0.00756
    p['rF'] = 0.35
    p['mF'] = 3.0
    p['IFxx'] = 0.1405
    p['IFyy'] = 0.28

    return p

def benchmark_to_moore(benchmarkParameters):
    """Returns the parameters for the Whipple model as derived by Jason K.
    Moore.

    Parameters
    ----------
    benchmarkParameters : dictionary
        Contains the set of parameters for the Whipple bicycle model as
        presented in Meijaard2007.

    Returns
    -------
    mooreParameters : dictionary
        The parameter set for the Moore derivation of the whipple bicycle model
        as presented in Moore2012.

    """
    bP = benchmarkParameters
    mP = {}

    # geometry
    mP['rf'] = bP['rF']
    mP['rr'] = bP['rR']
    mP['d1'] =  cos(bP['lam']) * (bP['c'] + bP['w'] - bP['rR'] * tan(bP['lam']))
    mP['d3'] = -cos(bP['lam']) * (bP['c'] - bP['rF'] * tan(bP['lam']))
    mP['d2'] = -cos(bP['lam'])*(bP['rF']-bP['rR']-bP['w']*tan(bP['lam']))

    # mass center locations
    # bicycle frame
    mP['l1'] = (bP['xB'] * cos(bP['lam']) - bP['zB'] * sin(bP['lam']) -
        bP['rR'] * sin(bP['lam']))
    mP['l2'] = (bP['xB'] * sin(bP['lam']) + bP['zB'] * cos(bP['lam']) +
        bP['rR'] * cos(bP['lam']))
    mP['l4'] = ((bP['zH'] + bP['rF']) * cos(bP['lam']) + (bP['xH'] - bP['w'])
        * sin(bP['lam']))
    mP['l3'] = ((bP['xH'] - bP['w'] - mP['l4'] * sin(bP['lam'])) /
        cos(bP['lam']))

    # masses
    mP['mc'] =  bP['mB']
    mP['md'] =  bP['mR']
    mP['me'] =  bP['mH']
    mP['mf'] =  bP['mF']

    # inertia
    # rear wheel inertia
    mP['id11']  =  bP['IRxx']
    mP['id22']  =  bP['IRyy']
    mP['id33']  =  bP['IRxx']

    # front wheel inertia
    mP['if11']  =  bP['IFxx']
    mP['if22']  =  bP['IFyy']
    mP['if33']  =  bP['IFxx']

    # lambda rotation matrix
    R = np.matrix([[cos(bP['lam']), 0, -sin(bP['lam'])],
            [0, 1,  0],          
            [sin(bP['lam']), 0,  cos(bP['lam'])]])

    # rotate the benchmark bicycle back frame inertia through the 
    #steer axis tilt, lambda
    IB =  np.matrix([[bP['IBxx'], 0., bP['IBxz']],
                     [0., bP['IByy'], 0.],
                     [bP['IBxz'], 0., bP['IBzz']]])
    IBrot =  R * IB * R.T

    # bicycle frame inertia
    mP['ic11'] =  IBrot[0, 0]
    mP['ic12'] =  IBrot[0, 1]
    mP['ic22'] =  IBrot[1, 1]
    mP['ic23'] =  IBrot[1, 2]
    mP['ic31'] =  IBrot[2, 0]
    mP['ic33'] =  IBrot[2, 2]


    # rotate the benchmark bicycle fork inertia through the steer axis tilt,
    # lambda
    IH =  np.matrix([[bP['IHxx'], 0., bP['IHxz']],
                     [0., bP['IHyy'], 0.],
                     [bP['IHxz'], 0., bP['IHzz']]])
    IHrot =  R * IH * R.T

    # fork/handlebar inertia
    mP['ie11'] =  IHrot[0, 0]
    mP['ie12'] =  IHrot[0, 1]
    mP['ie22'] =  IHrot[1, 1]
    mP['ie23'] =  IHrot[1, 2]
    mP['ie31'] =  IHrot[2, 0]
    mP['ie33'] =  IHrot[2, 2]

    # gravity
    mP['g'] = bP['g']

    return mP
