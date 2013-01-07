"""
Geometry description of steady turning:
1, turning radius of two wheels: Rr, Rf; (in model.py)
2, relative position of front wheel contact point fn to rear wheel dn:
   fn_dn_x, fn_dn_y; (in model.py)
3, the position of total mass center, tmc, of the bicycle relative to dn:
   tmc_dn_x, tmc_dn_y, Rc(turning radius of tmc);
4, alpha and beta angle according to the plot of steady turning for calculating
   lateral contact forces:
   Fx_r, Fy_r for rear wheel contact forces; Fx_f, Fy_f for front wheel.
   alpha: angle between vector of a_2 and lateral_v;
   beta: angle between vector of a_2 and centripetal direction.

All variables and other assumed names are shown in the plot.

FINALLY, ALPHA AND BETA ARE WHAT WE WANT FOR CALCULATION OF CONTACT FORCES.
"""

import numpy as np

def total_com(coordinates, masses):
    '''Returns the center of mass of a group of objects if the indivdual 
    position of centers of mass and mass is provided.

    coordinates : ndarray, shape(3,n)
        The rows are the x, y and z coordinates, respectively and the columns
        are for each object.
    masses : ndarray, shape(3,)
        An array of the masses of multiple objects, the order should correspond
        to the columns of coordinates.

    Returns
    -------
    mT : float
        Total mass of the objects.
    cT : ndarray, shape(3,)
        The x, y, and z coordinates of the total center of mass.

    '''
    products = masses * coordinates
    mT = np.sum(masses)
    cT = np.sum(products, axis=1) / mT
    return mT, cT
