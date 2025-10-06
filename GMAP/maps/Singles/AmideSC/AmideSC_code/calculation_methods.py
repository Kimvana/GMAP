
# 3rd party lib imports
from numba import njit
import numpy as np

# GMAP imports
import GMAP.src.tools.math_functions as GM_mf


def calc_dipole_Torii(map_, system, osc):
    """Calculate the Torii dipole moment for a given oscillator

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.SingleMap`
        The map instance which this function will belong to.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The class containing all the information on the system of the
        MD trajectory.
    osc : :class:`~GMAP.src.tools.system_reader.Oscillator`
        The specific oscillator for which the transformation is required.

    Returns
    -------
    r_vec : `np.ndarray`
        A length-3 vector containing the direction of the dipole moment.
        The vector must be normalized.
    r_pos : `np.ndarray`
        A length-3 vector containing the position of the dipole moment.
        The vector must lie within the simulation box.
    """

    pos_C_box = osc.positions_box[0]
    pos_O_box = osc.positions_box[1]
    pos_N_box = osc.positions_box[3]
    COvec = GM_mf.PBC_boxdiff_triclin(pos_O_box, pos_C_box, system.boxvects)
    COvec /= GM_mf.vec3_len(COvec)
    CNvec = GM_mf.PBC_boxdiff_triclin(pos_N_box, pos_C_box, system.boxvects)
    CNvec /= GM_mf.vec3_len(CNvec)

    # position of dipole vector
    r_pos = osc.positions[0] + 0.665*COvec + 0.258*CNvec

    # dipole moment vector itself
    r_vec = dipole_Torii(
        COvec, CNvec, map_.core.dipole_Torii_angle, map_.core.dipole_gas_phase)

    return r_vec, r_pos


@njit
def dipole_Torii(COvec, CNvec, itheta, magnitude):
    """Calculates the dipole moment using the Torii method.

    Method taken from AIM.

    Parameters
    ----------
    COvec : `np.ndarray`
        The vector pointing from the carbon to the oxygen atom
    CNvec : `np.ndarray`
        The vector pointing from the carbon to the nitrogen atom
    magnitude : `np.float32`
        The length that the vector is supposed to have.

    Returns
    -------
    mi : `np.ndarray`
        A length-3 vector containing the direction of the dipole moment.
        The vector must be normalized.
    """

    # itheta = 1/0.17632698  ## 1/tan(10 degrees expressed in radians)
    # itheta = 5.6712818196

    dri = 0.665*COvec + 0.258*CNvec

    dridri = GM_mf.dotprod(dri, dri)
    COvecdri = GM_mf.dotprod(COvec, dri)
    mi = dri - (COvecdri + np.sqrt(dridri - COvecdri*COvecdri)*itheta)*COvec

    mi /= GM_mf.vec3_len(mi)
    mi *= magnitude

    return mi


def get_position(map_, system, osc):
    """Determine the position for a given oscillator

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.SingleMap`
        The map instance which this function will belong to.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The class containing all the information on the system of the
        MD trajectory.
    osc : :class:`~GMAP.src.tools.system_reader.Oscillator`
        The specific oscillator for which the transformation is required.

    Returns
    -------
    position : `np.ndarray`
        A length-3 vector containing the position of the oscillator.
    """

    # position is correct, but shifted to wrong box
    pos = map_.code.GM_get_position_DMF(map_, system, osc)

    # now, shift it to the correct box
    boxpos = pos @ system.boxvects_inv
    boxpos -= np.floor(boxpos)

    return boxpos @ system.boxvects
