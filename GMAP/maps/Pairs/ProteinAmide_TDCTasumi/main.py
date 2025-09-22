
# 3rd party lib imports
from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.constants as GM_con
import GMAP.src.tools.math_functions as GM_mf


def GM_prep_coupling(map_, system, oscixlist, osclist):
    """Any preparation needed for calculating couplings this frame.

    In this case:
    When calculating the couplings, the dipole moment and position of
    each oscillator is needed for each combination of oscillators.
    Instead of calculating it again for each coupling, we do it once per
    oscillator here, so it can be read/reused often.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    oscixlist : list of int
        The oscillator indices of all oscillators that are treated by
        this map. Some might be only in a single pair, others in many.
    osclist : list of :class:`~GMAP.src.tools.system_reader.Oscillator`
        All oscillators treated by this map.
    """

    itheta = map_.Core.dipole_Torii_angle
    magnitude = map_.RunPars.Torii_dipole_magnitude

    for oscix, osc in zip(oscixlist, osclist):
        COvec = GM_mf.PBC_boxdiff_triclin(
            osc.positions_box[1], osc.positions_box[0], system.boxvects)
        COvec /= GM_mf.vec3_len(COvec)
        CNvec = GM_mf.PBC_boxdiff_triclin(
            osc.positions_box[3], osc.positions_box[0], system.boxvects)
        CNvec /= GM_mf.vec3_len(CNvec)
        dri = 0.665 * COvec + 0.258 * CNvec

        map_.dipole_pos_arr[oscix] = osc.positions[0] + dri

        COvecDri = GM_mf.dotprod(COvec, dri)
        dip_vec = dri - COvec * (COvecDri + itheta * np.sqrt(
            GM_mf.dotprod(dri, dri) - COvecDri * COvecDri))
        map_.dipole_vec_arr[oscix] = (
            dip_vec / GM_mf.vec3_len(dip_vec) * magnitude)

    map_.dipole_pos_arr = map_.dipole_pos_arr @ system.boxvects_inv


def GM_calc_coupling(map_, system, hamiltonian):
    """Calculate all the couplings that should be determined by this map

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    hamiltonian : `np.ndarray`
        The hamiltonian of the full system. Consists of float32, has a
        column and a row for each oscillator.
    """

    # We'll use this value a lot. It is the same as GM_con.e2i4pieps_angcm,
    # but now in units of cm^-1 * ang^3 Deb^-2
    i4pieps = np.float32(
        GM_con.i4pieps * GM_con.Debye**2 * GM_con.J2cm / GM_con.angstrom**3)
    for pair in map_.allpairs:
        J = calc_coupling(
            *pair, map_.dipole_pos_arr, map_.dipole_vec_arr, system.boxvects,
            i4pieps)
        hamiltonian[pair[0], pair[1]] = J
        hamiltonian[pair[1], pair[0]] = J


@njit
# same as DipDip map, as this is just DipDip map with Torii dipoles
def calc_coupling(oscix1, oscix2, pos_arr, vec_arr, boxvects, i4pieps):
    """Calculates the coupling value for the spcific provided pair.

    njit'ted for extra speed.

    Parameters
    ----------
    oscix1, oscix2 : int
        The oscillator index of each of the oscillators in this pair
        that should be calculated.
    pos_arr : `np.ndarray`
        The positions of all oscillators. All, not just those treated
        by this map. Although those not treated are not read, so
        can be set to 0.
    vec_arr : `np.ndarray`
        The dipole moment of all oscillators. All, not just those
        treated by this map. Although those not treated are not read, so
        can be set to 0.
    boxvects : `np.ndarray`
        The box vectors of the MD system.

    Returns
    -------
    J : float
        The coupling found for the provided pair.
    """

    # Used constants:
    # Cm = (1/3.33564) * 10^30 D  (Coulomb meter in Debye)
    # m = 10^10 ang (meter in angstrom)
    # J = 1/hc = (1/1.98644586) * 10^25 1/m
    # => J = 5.03411656 * 10^22 1/cm (joule in wavenumbers)
    # eps_0 = 8.8541878128 F/m = 8.8541878128 C^2/Jm (coulomb squared per
    # joule meter)

    # derived value:
    # 4piEinv = 1/(4 * pi * eps_0) Jm/C^2
    # Gives 5034.11656 cm^-1 * ang^3 Deb^-2

    # the positions array is in box-coordinates -> easy subtraction, then
    # move back into cartesian
    d = GM_mf.PBC_back2box(pos_arr[oscix1, :] - pos_arr[oscix2, :], boxvects)
    ir2 = 1/GM_mf.dotprod(d, d)
    ir = np.sqrt(ir2)
    ir3 = ir*ir2
    ir5 = ir3*ir2

    return i4pieps * (
        GM_mf.dotprod(vec_arr[oscix1], vec_arr[oscix2]) * ir3
        - 3.0 * GM_mf.dotprod(vec_arr[oscix1], d)
        * GM_mf.dotprod(vec_arr[oscix2], d) * ir5)


def GM_pre_run(map_, system):
    """Initialize the data structure for saving dipole moments and
    positions.

    This is used so save prep_calc's preparation.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    """

    map_.dipole_vec_arr = np.zeros((system.nosc, 3), dtype="float32")
    map_.dipole_pos_arr = np.zeros((system.nosc, 3), dtype="float32")
    map_.Core.dipole_Torii_angle = np.float32(
        1 / np.tan(GM_con.deg2rad * map_.RunPars.Torii_dipole_angle))
