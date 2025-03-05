import numpy as np
import GMAP.src.tools.MathFunctions as GM_MF

"""
Overwrite the Centre of Mass with the position of the Mg atom in determining
which bChl-c molecules are within the radius of consideration.

This feature might not be needed beyond comparison with previous work.
"""

def GM_pre_frame(map_,system):
    mg_idx = 26
    for oscillator in system.oscillators_ordered["BChl-c"]:
        resnum = system.resnums[oscillator.used_atoms[mg_idx]]
        system.residues.CoM[resnum] = oscillator.positions_box[mg_idx]
