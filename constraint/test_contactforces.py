"""
test_contactforces module:
1, test contact forces in a given steady turning condition.
   one is derived from geometry of steady turning,
   the other is from KaneMethod in mech.
2, Actually, comparison should happen to between the manual derivation and 
   KaneMethod class, along with same defined variables.
"""

from contactforces_steadyturning import SteadyTurning
from bicycle import benchmark_parameters

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
    Fx_r, Fy_r, Fz_r, Fx_f, Fy_f, Fz_f = st._contactforcesSym
    u1 = st._equilibriumSym[0]
    yawrate = st._equilibrium_u[u1]

    alpha = arctan (float((fda1/(Rr - fda2))))
    beta = arctan (float((coma1/(Rr - coma2))))

    # Equations according to geometry.
    # Lateral forces
    tangential = sin(beta) * Fy_r - sin(alpha-beta) * Fy_f
    centrifugal = cos(alpha-beta) * Fy_f + cos(beta) * Fy_r - mT * yawrate**2 * R
    # Normal forces
    forcial = Fz_r + Fz_f - mT * benchmark_parameters()['g']
    torqueal = coma1 * Fz_r - (fda1 - coma1) * Fz_f

    # Calculate the contact forces in geometric way;
    # Calculate the contact forces from the model.
    contact_forces_geo = solve([tangential, centrifugal, forcial, torqueal], 
                               [Fy_r, Fz_r, Fy_f, Fz_f])
    Fy_r_geo = contact_forces_geo[Fy_r]
    Fz_r_geo = contact_forces_geo[Fz_r]
    Fy_f_geo = contact_forces_geo[Fy_f]
    Fz_f_geo = contact_forces_geo[Fz_f]

    contact_forces_model = st.contact_force()
    Fx_r_mod = contact_forces_model[Fx_r]
    Fy_r_mod = contact_forces_model[Fy_r]
    Fz_r_mod = contact_forces_model[Fz_r]
    Fx_f_mod = contact_forces_model[Fx_f]
    Fy_f_mod = contact_forces_model[Fy_f]
    Fz_f_mod = contact_forces_model[Fz_f]

    contact_forces_dict = {'geo [rearxyz, frontxyz]': [(Fy_r_geo, Fz_r_geo), 
                                                       (Fy_f_geo, Fz_f_geo)],
                           'mod [rearxyz, frontxyz]': [(Fx_r_mod, Fy_r_mod, Fz_r_mod),
                                                       (Fx_f_mod, Fy_f_mod, Fz_f_mod)]}

    # compare
    #assert (Fy_r_geo == Fy_r_mod)
    #assert (Fy_f_geo == Fy_f_mod)

    return contact_forces_dict
