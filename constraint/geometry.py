"""
Geometry description of steady turning:
1, turning radius of two wheels: Rr, Rf;
2, relative position of front wheel contact point fn to rear wheel dn:
   fn_dn_x, fn_dn_y;
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

