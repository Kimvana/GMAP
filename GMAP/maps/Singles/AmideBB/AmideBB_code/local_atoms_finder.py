

def find_local_atoms(map_, system, oscillator_list):
    """Finds the local atoms for all oscillators in this map.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    oscillator_list : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        All oscillators belonging to a single struct of this map.
    """

    if map_.RunPars.frequency_map_choice == "Tokmakoff":
        get_Tokmakoff_locals(oscillator_list)

        # bring updates to c-friendly versions of the local_atoms list
        for oscillator in oscillator_list:
            oscillator.process_atschoice("local_atoms")
        return  # nothing else should happen for Tokmakoff

    add_proline_HDs(system, oscillator_list)

    if map_.RunPars.consider_nearest_neighbours:
        for osc in oscillator_list:
            if osc.NtermNB is not None:
                osc.local_atoms.extend(get_Nterm_locals(system, osc.NtermNB))
            if osc.CtermNB is not None:
                osc.local_atoms.extend(get_Cterm_locals(osc.CtermNB))
            osc.local_atoms.extend(get_proline_atoms(system, osc))

    for osc in oscillator_list:
        osc.local_atoms.extend(get_CA_hydrogen_atoms(map_, system, osc))

    # bring updates to c-friendly versions of the local_atoms list
    for oscillator in oscillator_list:
        oscillator.process_atschoice("local_atoms")


def get_Tokmakoff_locals(oscillator_list):
    """Get all atoms that are considered local for the Tokmakoff map.

    Parameters
    ----------
    oscillator_list : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        All oscillators belonging to a single struct of this map.
    """

    # For Tokmakoff, the locals == used_atoms (as in core.txt). The only
    # exception to this is the prepros, they shouldn't contain the CD atom.
    for osc in oscillator_list:
        if osc.resnames[1] == "PRO":  # prepro!
            osc.local_atoms = osc.used_atoms[:4] + [osc.used_atoms[5]]


def add_proline_HDs(system, oscillator_list):
    """Get all HD atoms for all prolines in the system.

    Parameters
    ----------
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    oscillator_list : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        All oscillators belonging to a single struct of this map.
    """

    sysres = system.residues
    for osc in oscillator_list:
        if osc.resnames[1] == "PRO":
            resix = osc.resnums[1]
            firstix = sysres.first_ix[resix]
            lastix = sysres.last_ix[resix]

            # atoms in this residue
            for atix, atname in zip(
                system.atnums[firstix:lastix + 1],
                system.atnames[firstix:lastix + 1]
            ):
                if atname in ("HD1", "HD2"):
                    osc.local_atoms.append(atix)


def get_Nterm_locals(system, osc):
    """Get all local atoms on the Nterm side of the oscillator

    Parameters
    ----------
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The oscillator for which these local atoms need to be found
    """

    # for the opls forcefield, the non-CA atoms of the N term neighbour should
    # be added to locals. Otherwise, all atoms of the N term neighbour should
    # be added to locals
    if system.types[0][:4] == "opls":
        return [osc.used_atoms[ix] for ix in [0, 1, 3, 4]]
    else:
        return osc.used_atoms


def get_Cterm_locals(osc):
    """Get all local atoms on the Cterm side of the oscillator

    Parameters
    ----------
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The oscillator for which these local atoms need to be found
    """

    # ETA residue is different (not yet implemented in detection)
    if osc.resnames[1] == "ETA":
        local_ix = osc.used_atoms[0:5]
    else:
        local_ix = osc.used_atoms
    return local_ix


def get_proline_atoms(system, osc):
    """Add all atoms in neighbouring prolines to the locals

    Parameters
    ----------
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The oscillator for which these local atoms need to be found
    """

    # any atoms in proline residues should be taken to be local, too. This
    # could be either of this oscillators residues, or the one after.
    resnums = list(osc.resnums)
    if osc.CtermNB is not None:
        resnums.append(osc.CtermNB.resnums[1])
    proline_atoms = []
    for resnum in resnums:
        if system.residues.resnames[resnum] == "PRO":
            proline_atoms.extend(range(
                system.residues.first_ix[resnum],
                system.residues.last_ix[resnum] + 1
            ))

    # But if the one after is a proline, its C=O group shouldn't be added.
    if osc.CtermNB is not None and osc.CtermNB.resnames[1] == "PRO":
        resix = osc.CtermNB.resnums[1]
        firstix = system.residues.first_ix[resix]
        lastix = system.residues.last_ix[resix]

        # atoms in this residue
        for atix, atname in zip(
            system.atnums[firstix:lastix + 1],
            system.atnames[firstix:lastix + 1]
        ):
            if atname in ("C", "O", "OC1", "OC2"):
                proline_atoms.remove(atix)

    return proline_atoms


def get_CA_hydrogen_atoms(map_, system, osc):
    """Get all hydrogen atoms on all CA atoms already in locals

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The oscillator for which these local atoms need to be found
    """

    CA_hydrogen_atoms = []
    resnums = list(osc.resnums)
    if map_.RunPars.consider_nearest_neighbours:
        if osc.CtermNB is not None:
            resnums.append(osc.CtermNB.resnums[1])
        if osc.NtermNB is not None and system.types[0][:4] != "opls":
            resnums.append(osc.NtermNB.resnums[0])

    for resnum in resnums:
        firstix = system.residues.first_ix[resnum]
        lastix = system.residues.last_ix[resnum]

        # atoms in this residue
        for atix, atname in zip(
            system.atnums[firstix:lastix + 1],
            system.atnames[firstix:lastix + 1]
        ):
            if atname in ("HA", "HA1", "HA2", "HA3"):
                CA_hydrogen_atoms.append(atix)

    return CA_hydrogen_atoms
