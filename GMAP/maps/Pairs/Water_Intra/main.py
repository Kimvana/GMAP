
# 3rd party lib imports
# from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.MathFunctions as GM_MF
from GMAP.src.tools.PrintTools import devprint as dpr


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

    return

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
    x1=0.1934-1.75e-5*osc1.freq
    x2=0.1934-1.75e-5*osc2.freq
    p1=1.611+5.893e-4*osc1.freq
    p2=1.611+5.893e-4*osc2.freq
    sumE=osc1.VEGout[0,1]+osc2.VEGout[0,1]
    J=(-1789+23852*sumE)*x1*x2-1.966*p1*p2

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
