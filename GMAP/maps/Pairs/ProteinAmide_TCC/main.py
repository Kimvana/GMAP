
# 3rd party lib imports
# from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.MathFunctions as GM_MF
# from GMAP.src.tools.PrintTools import devprint as dpr

# own module imports
import ProteinAmide_TCC_code.TCCclib as MC_TC


def GM_prep_coupling(map_, system, oscixlist, osclist):
    """Any preparation needed for calculating couplings this frame.

    In this case:
    When calculating the couplings, this property 'v' of each oscillator
    is needed for each combination of oscillators. Instead of
    calculating it again for each coupling, we do it once per oscillator
    here, so it can be read/reused often.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    oscixlist : list of int
        The oscillator indices of all oscillators that are treated by
        this map. Some might be only in a single pair, others in many.
    osclist : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        All oscillators treated by this map.
    """

    for oscix, osc in zip(oscixlist, osclist):
        COvec = GM_MF.PBC_boxdiff_triclin(
            osc.positions_box[1], osc.positions_box[0], system.boxvects)
        COvec /= GM_MF.vec3_len(COvec)

        CNvec = GM_MF.PBC_boxdiff_triclin(
            osc.positions_box[3], osc.positions_box[0], system.boxvects)
        CNvec = GM_MF.project(COvec, CNvec)
        CNvec /= GM_MF.vec3_len(CNvec)
        z = GM_MF.crossprod(COvec, CNvec)
        z /= GM_MF.vec3_len(z)

        if osc.Map.name == "AmideSC" or osc.resnames[1] != "PRO":
            alpha = map_.alpha_gen
            v = map_.v_gen
        else:
            alpha = map_.alpha_pro
            v = map_.v_pro

        rotmat = np.array([COvec, CNvec, z])
        map_.map_tcc_v[oscix] = np.dot(v, rotmat) * alpha


def GM_calc_coupling(map_, system, hamiltonian):
    """Calculate all the couplings that should be determined by this map

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    hamiltonian : `np.ndarray`
        The hamiltonian of the full system. Consists of float32, has a
        column and a row for each oscillator.
    """

    for pair in map_.allpairs:
        J = calc_coupling(*pair, map_, system)
        hamiltonian[pair[0], pair[1]] = J
        hamiltonian[pair[1], pair[0]] = J


def calc_coupling(oscix1, oscix2, map_, system):
    """Calculates the coupling value for the spcific provided pair.

    Parameters
    ----------
    oscix1, oscix2 : int
        The oscillator index of each of the oscillators in this pair
        that should be calculated.
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everyting the program currently knows
        about the MD system.

    Returns
    -------
    J : float
        The coupling found for the provided pair.
    """

    osc1 = system.oscillators[oscix1]
    osc2 = system.oscillators[oscix2]

    if osc1.Map.name == "AmideSC" or osc1.resnames[1] != "PRO":
        q1 = map_.q_gen
        dq1 = map_.dq_gen
    else:
        q1 = map_.q_pro
        dq1 = map_.dq_pro

    if osc2.Map.name == "AmideSC" or osc2.resnames[1] != "PRO":
        q2 = map_.q_gen
        dq2 = map_.dq_gen
    else:
        q2 = map_.q_pro
        dq2 = map_.dq_pro

    diff = osc2.positions_box[None, :, :] - osc1.positions_box[:, None, :]
    diff -= np.floor(diff + 0.5)
    diff = diff @ system.boxvects
    r2 = np.sum(diff * diff, axis=2)

    r2[r2 < 0.01] = 1
    ir2 = 1 / r2
    ir = np.sqrt(ir2)
    ir3 = ir * ir2
    ir5 = ir3 * ir2

    v1 = map_.map_tcc_v[oscix1]
    v2 = map_.map_tcc_v[oscix2]
    # dpr(oscix1, oscix2, "\n", v1, "\n", v2)
    # dpr("XXXXXXXXXX", np.sum(v2[None, :, :] * diff, axis=2))

    J = ir * dq1[:, None] * dq2[None, :]
    J -= ir3 * (
        dq1[:, None] * q2[None, :] * np.sum(v2[None, :, :] * diff, axis=2)
        - q1[:, None] * dq2[None, :] * np.sum(v1[:, None, :] * diff, axis=2)
        - np.sum(v1[:, None, :] * v2[None, :, :], axis=2)
        * q1[:, None] * q2[None, :])
    J -= (
        3 * ir5 * q1[:, None] * q2[None, :]
        * np.sum(v2[None, :, :] * diff, axis=2)
        * np.sum(v1[:, None, :] * diff, axis=2))
    J = np.sum(J)
    J *= map_.fourPiEps

    return J


def GM_pre_run(map_, system):
    """Initialize the data structure for saving v.

    This is used so save prep_calc's preparation.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    """

    map_.map_tcc_v = np.zeros((system.nosc, 6, 3), dtype="float32")


def GM_post_init(map_, system):
    """Do some final initializations that need to happen before the
    calculation starts.

    Steps present:
    - Extract all map parameters

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    """

    # read in all parameters from the constants file
    with open(map_.directory / "constants.txt") as fhand:
        # get_next_line just gets the next not-empty line from the file
        map_.fourPiEps = float(get_next_line(fhand).split()[0])
        map_.alpha_gen = float(get_next_line(fhand).split()[0])
        map_.alpha_pro = float(get_next_line(fhand).split()[0])

        map_.q_gen = np.array(
            [float(item) for item in get_next_line(fhand).split()])
        map_.q_pro = np.array(
            [float(item) for item in get_next_line(fhand).split()])
        map_.dq_gen = np.array(
            [float(item) for item in get_next_line(fhand).split()])
        map_.dq_pro = np.array(
            [float(item) for item in get_next_line(fhand).split()])

        map_.v_gen = np.array([
            [float(item) for item in get_next_line(fhand).split()]
            for _ in range(6)])
        map_.v_pro = np.array([
            [float(item) for item in get_next_line(fhand).split()]
            for _ in range(6)])

    # We'd like to use the c-library for this map (so it is considerably
    # faster)
    MC_TC.init_map_for_clib(map_, system)


def get_next_line(fhand):
    """Returns the next non-empty line from the provided file handle.

    This method ignores comments in the file - so a line with only
    a comment is considered an empty line.

    Parameters
    ----------
    fhand : `_io.TextIOWrapper`
        The file (handle) from which the next line is desired

    Returns
    -------
    line : str
        The next non-empty line from the file. Comments are ignored when
        parsing the file.
    """

    line = ""
    while len(line) == 0:
        line = fhand.readline()
        line = line.split("#")[0].strip()
    return line
