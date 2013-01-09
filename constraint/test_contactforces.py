"""
test_contactforces module:
1, test contact forces in a given steady turning condition.
   here, calculate the contact forces directly from geometry.
"""

from contactforces_steadyturning import SteadyTurning

from numpy import (sin, cos, arctan)
from sympy import solve
import pdb

# Import SteadyTurning class and build geometry of a given steady turning.
# Import alpha and beta angles for calculation of contact forces.
# Contact forces also have normal forces on each wheel.
def contact_forces_test(lean, steer):
    """Return contact forces."""

    # Import SteadyTurning class
    st = SteadyTurning(lean, steer)

    Rr = st._turningRadiusRearGeo
    R = st._turningRadiusCom
    mT = st._totalMass
    fda1 = st._fn_dn_A1
    fda2 = st._fn_dn_A2
    coma1 = st._totalComA123[0]
    coma2 = st._totalComA123[1]
    Fy_r = st._contactforcesSym[1]
    Fy_f = st._contactforcesSym[3]
    u1 = st._equilibriumSym[0]
    yawrate = st._equilibrium_u[u1]

    alpha = arctan (float((fda1/(Rr - fda2))))
    beta = arctan (float((coma1/(Rr - coma2))))

    tangential = sin(beta) * Fy_r - sin(alpha-beta) * Fy_f
    centrifugal = cos(alpha-beta) * Fy_f + cos(beta) * Fy_r - mT * yawrate**2 * R

    # Calculate the contact forces in geometric way;
    # Calculate the contact forces from the model.
    contact_forces_geo = solve([tangential, centrifugal], [Fy_r, Fy_f])
    Fy_r_geo = contact_forces_geo[Fy_r]
    Fy_f_geo = contact_forces_geo[Fy_f]

    contact_forces_model = st.contact_force()
    Fy_r_mod = contact_forces_model[Fy_r]
    Fy_f_mod = contact_forces_model[Fy_f]

    contact_forces_dict = {'geo (rearxyz, frontxyz)': (Fy_r_geo, Fy_f_geo),
                           'mod (rearxyz, frontxyz)': (Fy_r_mod, Fy_f_mod)}
    # compare
    #assert (Fy_r_geo == Fy_r_mod)
    #assert (Fy_f_geo == Fy_f_mod)

    return contact_forces_dict
