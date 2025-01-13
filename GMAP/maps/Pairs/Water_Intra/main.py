
# gmap imports
import GMAP.src.tools.constants as GM_Con


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

    # These are the couplings according to the OH intramolecular water map
    # We follow the equations of Skinner to get gecognizable parameters
    x1 = 0.1934 - 1.75e-5 * osc1.freq
    x2 = 0.1934 - 1.75e-5 * osc2.freq
    p1 = 1.611 + 5.893e-4 * osc1.freq
    p2 = 1.611 + 5.893e-4 * osc2.freq
    # Find the sum of the fields in atomic units
    sumE = (osc1.VEGout[0, 1] + osc2.VEGout[0, 1]) * GM_Con.bohr2ang ** 2
    J = (-1789 + 23852 * sumE) * x1 * x2 - 1.966 * p1 * p2

    return J
