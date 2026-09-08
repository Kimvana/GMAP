"""
This file is NOT being used by GMAP. If it should be used, this file
should be renamed to "main.py" (current name is main_MG_position.py).

The purpose of this file is to replace all centres of mass of all BChlC
residues with the position of their magnesium atom. This means that
when determining the influencers for electrostatics calculations, some
residues might fall in/out of the sphere of influence. The sphere itself
is NOT affected by this.
"""


def GM_pre_frame(map_, system):
    """Overwrite the Centre of Mass with the position of the Mg atom in
    determining which bChl-c molecules are within the radius of
    consideration.

    This feature might not be needed beyond comparison with previous
    work.
    """

    mg_idx = 26
    for oscillator in system.oscillators_ordered["BChl_c"]:
        resnum = system.resnums[oscillator.used_atoms[mg_idx]]
        system.residues.CoM[resnum] = oscillator.positions_box[mg_idx]
