

def GM_adjust_oscillators(map_, system, oscillator_list):
    """Makes the necessary changes to the list of oscillators.

    The program finds all oscillators mathing the instructions from
    core.txt. However, there is no way for the program to avoid double
    counting symmetrical groups (like the cystbridge mockup example).
    If a map knows its group is symmetrical, this function can be
    designed to only return half of the inputs.

    Another possible use is for the code of the map to get to know its
    oscillators. When all oscillators are passed through this function,
    the (global) atom number of the first atom of this group (for
    example) can be linked to a specific property the group might need
    to know. This might be useful if a map needs to cover two very
    similar oscillators.

    .. note::
        This function is called separately for each struct that the map
        defines. So take into account that the function could be called
        multiple times within a single simulation!

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    oscillator_list : list of :class:`~GMAP.src.tools.system_reader.Oscillator`
        All oscillators belonging to a single struct of this map.

    Returns
    -------
    oscillator_list : list of :class:`~GMAP.src.tools.system_reader.Oscillator`
        All oscillators belonging to a single struct of this map.
    """

    for oscillator in oscillator_list:
        print(oscillator)
        print(oscillator.used_atoms)

    return oscillator_list
