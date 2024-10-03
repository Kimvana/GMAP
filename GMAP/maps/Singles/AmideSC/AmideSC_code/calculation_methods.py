
# 3rd party lib imports
from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.MathFunctions as GM_MF


# def calc_dipole_Jansen(map_, system, osc):
#     map_.Core.dipole_gas_phase_array   # gas phase
#     map_.Core.dipole_data_array   # arrays.

#     r_pos = osc.positions[0] + 0.868 * osc.rotation_matrix[0]

#     r_vec = map_.Core.dipole_gas_phase_array + np.sum(
#         np.multiply(osc.VEGout[None, :, :], map_.Core.dipole_data_array),
#         axis=(1, 2)
#     )

#     return r_vec, r_pos


def calc_dipole_Torii(map_, system, osc):
    rotmat = osc.rotation_matrix
    COvec = rotmat[0]
    CNvec = rotmat[1]

    # position of dipole vector
    r_pos = osc.positions[0] + 0.665*COvec + 0.258*CNvec

    # dipole moment vector itself
    r_vec = dipole_Torii(COvec, CNvec, map_.Core.dipole_gas_phase)

    return r_vec, r_pos


@njit
def dipole_Torii(COvec, CNvec, magnitude):
    """Calculates the dipole moment using the Torii method.

    Method taken from AIM.
    """

    # itheta = 1/0.17632698  ## 1/tan(10 degrees expressed in radians)
    itheta = 5.6712818196

    dri = 0.665*COvec + 0.258*CNvec

    dridri = GM_MF.dotprod(dri, dri)
    COvecdri = GM_MF.dotprod(COvec, dri)
    mi = dri - (COvecdri + np.sqrt(dridri - COvecdri*COvecdri)*itheta)*COvec

    mi /= GM_MF.vec3_len(mi)
    mi *= magnitude

    return mi
