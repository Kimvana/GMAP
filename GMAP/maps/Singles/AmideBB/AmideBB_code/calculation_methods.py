
# 3rd party lib imports
from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.MathFunctions as GM_MF
import GMAP.src.tools.PrintTools as GM_PT


def calc_dipole_Torii(map_, system, osc):
    """Calculate the Torii dipole moment for a given oscillator

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.SingleMap`
        The map instance which this function will belong to.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
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
    itheta = 5.6712818196

    dri = 0.665*COvec + 0.258*CNvec

    dridri = GM_MF.dotprod(dri, dri)
    COvecdri = GM_MF.dotprod(COvec, dri)
    mi = dri - (COvecdri + np.sqrt(dridri - COvecdri*COvecdri)*itheta)*COvec

    mi /= GM_MF.vec3_len(mi)
    mi *= magnitude

    return mi


def calc_dipole_Jansen(map_, system, osc):
    """Calculate the Jansen dipole moment for a given oscillator

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.SingleMap`
        The map instance which this function will belong to.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
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

    _, r_pos = map_.code.GM_get_dipole_dir(map_, system, osc)
    if osc.resnames[1] == "PRO":
        gasdip = map_.Core.dipole_gas_phase_array_prepro
        diparr = map_.Core.dipole_data_array_prepro
    else:
        gasdip = map_.Core.dipole_gas_phase_array
        diparr = map_.Core.dipole_data_array
    print(osc.oscix, gasdip)
    xyz_local = gasdip + np.sum(
        np.multiply(osc.VEGout[None, :, :], diparr), axis=(1, 2))
    xyz_cartesian = np.dot(xyz_local, osc.rotation_matrix)
    return xyz_cartesian, r_pos


def get_position(map_, system, osc):
    """Determine the position for a given oscillator

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.SingleMap`
        The map instance which this function will belong to.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
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


def neighbor_influence(map_, system, osc):
    """Calculates how much the frequency should be shifted due to neighbors

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.SingleMap`
        The map instance which this function will belong to.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The specific oscillator for which the transformation is required.

    Returns
    -------
    delta : float
        The shift that should be added to the frequency
    """

    delta = 0
    if osc.NtermNB is not None:
        delta += map_.neighbormaps["NtermShift" + osc.Nterm_nnmap].get_delta(
            osc.NtermNB, osc, system)
    if osc.CtermNB is not None:
        delta += map_.neighbormaps["CtermShift" + osc.Cterm_nnmap].get_delta(
            osc, osc.CtermNB, system)
    return delta


def determine_maps(oscillator_list, map_, system):
    """
    For each oscillator, find out what map should be used to consider its
    N-term and C-term neighbor.

    The determined maps are saved as an attribute of the oscillators.

    Parameters
    ----------
    oscillator_list : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        All oscillators belonging to a single struct of this map.
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    """

    for osc in oscillator_list:
        # N term is first, C term is last
        if osc.NtermNB is not None:
            # supply both oscs.
            osc.Nterm_nnmap = determine_map(osc.NtermNB, osc, map_, system)
        if osc.CtermNB is not None:
            osc.Cterm_nnmap = determine_map(osc, osc.CtermNB, map_, system)


def determine_map(osc1, osc2, map_, system):
    """Determines which nearest-neighbor map should be used for this
    pair of oscillators

    Parameters
    ----------
    osc1, osc2 : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The specific oscillators for which the neighbour mapping is
        determined.
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)

    Returns
    -------
    mapname : str
        The name (without 'Nterm' or 'Cterm' prefix) of the map that
        fits the given pair.
    """

    # no (pre-)prolines! easy!
    if "PRO" not in osc2.resnames:
        return ""

    # for now, the pro-pro case is treated as if it is pro-gly. in the future,
    # a pro-pro map should be made!

    boxpos = osc1.positions_box
    COvec = GM_MF.PBC_boxdiff_triclin(boxpos[1], boxpos[0], system.boxvects)
    NHvec = GM_MF.PBC_boxdiff_triclin(boxpos[4], boxpos[3], system.boxvects)

    # (both pro-pro(should for now be treated as pro-gly) and pro-gly)
    if osc1.resnames[1] == "PRO":
        # In the original code, Pro-Pro is actually treated as Gly-Pro
        if (
            osc2.resnames[1] == "PRO"
            and map_.RunPars.legacy_mode == "AmideImaps"
        ):
            bondtype = "GP"
        else:
            bondtype = "PG"
    else:
        bondtype = "GP"

    if bondtype == "GP":
        if GM_MF.dotprod(COvec, NHvec) < 0:
            return "_transGly_transPro"
        else:
            return "_cisGly_transPro"

    LorD = DLcheck(osc1, osc2, map_, system)  # -1 for D, 1 for L, 0 for nodir
    if GM_MF.dotprod(COvec, NHvec) < 0:
        if LorD < 0:
            return "_transDPro_transGly"
        else:
            return "_transPro_transGly"
    elif LorD < 0:
        return "_cisDPro_transGly"
    else:
        return "_cisPro_transGly"


def DLcheck(osc1, osc2, map_, system):
    """Determine of an amino acid is of D or L chirality

    Parameters
    ----------
    osc1, osc2 : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The oscillators surrounding the amino acid in question
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)

    Returns
    -------
    chirality : int
        Whether the amino acid is in D or L configuration.
        -1 is returned for D (special), 1 is returned for L (default for
        most organisms), 0 is returned for groups where the CA atom
        is not chiral (like glycine for example)
    """

    # checks whether the amino acid between two oscillators is in L or D
    # configuration

    if osc2.resnames[0] in ("GLY", "FOR", "ETA", "GL2"):
        return 0

    atomCA = osc2.positions_box[2]
    atomC = osc2.positions_box[0]
    atomN = osc1.positions_box[3]

    first_at = system.residues.first_ix[osc2.resnums[0]]
    last_at = system.residues.last_ix[osc2.resnums[0]]
    for atix, atname in zip(
        system.atnums[first_at:last_at+1],
        system.atnames[first_at:last_at+1]
    ):
        if atname == "CB":
            atomCBix = atix
            break
    else:
        GM_PT.Printer.warning(
            "Warning! The residue between the following two oscillators "
            "does not have a CB atom, and thus its chirality cannot be "
            f"determined:\n{osc1}\n{osc2}\nPlease make sure you're applying "
            "the correct map to the correct system. If this map erraneously "
            "detects something it shouldn't, please contact the developers!",
            "map_AmideBB_4"
        )
        map_.success = False
        return 0
    shift = atomC
    atomCB = system.positions[atomCBix] @ system.boxvects_inv - shift
    atomCB -= np.floor(atomCB + 0.5) - shift

    # used ats order:    res0{C O CA} res1{N H CA} ({N CD CA} for prepro)
    # osc1 is first, osc2 is last
    CACvec = GM_MF.PBC_boxdiff_triclin(atomC, atomCA, system.boxvects)
    CANvec = GM_MF.PBC_boxdiff_triclin(atomN, atomCA, system.boxvects)
    CACBvec = GM_MF.PBC_boxdiff_triclin(atomCB, atomCA, system.boxvects)

    CxN = GM_MF.crossprod(CACvec, CANvec)
    if GM_MF.dotprod(CxN, CACBvec) > 0:
        return -1
    else:
        return 1
