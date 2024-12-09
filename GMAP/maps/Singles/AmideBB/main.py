
# 3rd party imports
import numpy as np

# GMAP imports
import GMAP.src.tools.constants as GM_con
import GMAP.src.tools.PrintTools as GM_PT
# from GMAP.src.tools.PrintTools import devprint as dpr

# own module imports
import AmideBB_code.calculation_methods as MC_CM
import AmideBB_code.local_atoms_finder as MC_LAF
import AmideBB_code.neighbor_manager as MC_NM
import AmideBB_code.parameter_changer as MC_PC


def GM_adjust_map_core_raw(map_):
    """Makes the necessary changes to the 'raw' input read from core.txt.

    Is expected to not return anything - return value is not caught.

    The core.txt file is stored in Map.rawcore. It has not yet been
    parsed, just loaded into a dictionary. In this dictionary, each
    keyword is its own dictionary key. Most keywords can only occur once
    in the file - those have a list of the 'words' on the line as
    their value. The parameters that are allowed to occur more than once
    have a list as value, in which other lists appear - one for each
    line.

    The core.txt file has to be changed because the parameters of this
    map allow to change between models, each of which has their own
    files.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    # adjusting the names of functional_group
    all_amino_acid_codes = [
        "ARG", "HIS", "LYS", "ASP", "GLU", "SER", "THR", "ASN", "GLN",
        "CYS", "GLY", "PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE",
        "TYR", "TRP"
    ]
    amino_acids_joined = ",".join(all_amino_acid_codes)

    oldentry = map_.rawcore["functional_group"]

    map_.rawcore["functional_group"] = [
        [
            word.replace("anyprot", amino_acids_joined)
            for word in struct
        ] for struct in oldentry
    ]

    MC_PC.adjust_map_core_raw(map_)


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
    Map : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    Syst : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    oscillator_list : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        All oscillators belonging to a single struct of this map.

    Returns
    -------
    oscillator_list : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        All oscillators belonging to a single struct of this map.
    """

    # Sort oscillators into the correct order
    oscillator_list = MC_PC.oscillator_sorter(map_, system, oscillator_list)

    # tell each oscillator what/who it's neighbors are.
    # N term is first, C term is last
    for oscillator in oscillator_list:
        oscillator.NtermNB = None
        oscillator.CtermNB = None
    # used ats order:    res0{C O CA} res1{N H CA} ({N CD CA} for prepro)
    for Nosc in oscillator_list:
        for Cosc in oscillator_list:
            if Nosc.used_atoms[5] == Cosc.used_atoms[2]:
                Nosc.CtermNB = Cosc
                Cosc.NtermNB = Nosc
                break  # Nosc can at most have a single Cterm neighbour

    return oscillator_list


# A place to do further initialization if a map requires it. Think of
# things like building further lookup tables, for instance.
# (for AmideBB - find neighbours!)
def GM_post_init(map_, system):
    """Do some final initializations that need to happen before the
    calculation starts.

    Checks include:
    - comparing RunPars of this map to that of AmideSC, if the latter is
      present and active
    - initializing the prepro properties/files
    - Assigning the correct functions based on the parameter choices
    - Finding and assigning the neighbours of each group
    - Identifying all atoms local to each oscillator.

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

    # verify that both dipole maps (if applicable) have the same map choices
    rps = map_.RunPars
    main_runpars = map_.RunPars.MainRunPars
    # is other map present?
    if "AmideSC" in main_runpars.requested_mapdict.keys():
        amSC_rps = main_runpars.requested_mapdict["AmideSC"].RunPars
        if (  # the maps do not match, and they're not allowed to mismatch.
            amSC_rps.frequency_map_choice != rps.frequency_map_choice
            and not rps.allow_map_mismatch
        ):
            GM_PT.Printer.warning(
                "Warning! The current calculation makes use of both the "
                "AmideBB and AmideSC maps, but they make use of different "
                "frequency maps. For most physical applications, this does "
                "not make sense. If you are absolutely sure that you want "
                "the two maps to have different methods, make sure to set "
                "AmideBB.allow_map_mismatch to true. When in doubt, "
                "consult the README.",
                "map_AmideBB_2", True
            )
        if (  # the maps do not match, and they're not allowed to mismatch.
            amSC_rps.dipole_map_choice != rps.dipole_map_choice
            and not rps.allow_map_mismatch
        ):
            GM_PT.Printer.warning(
                "Warning! The current calculation makes use of both the "
                "AmideBB and AmideSC maps, but they make use of different "
                "dipole maps. For most physical applications, this does "
                "not make sense. If you are absolutely sure that you want "
                "the two maps to have different methods, make sure to set "
                "AmideBB.allow_map_mismatch to true. When in doubt, "
                "consult the README.",
                "map_AmideBB_2", True
            )
        if (  # the maps do not match, and they're not allowed to mismatch.
            amSC_rps.legacy_mode != rps.legacy_mode
            and not rps.allow_map_mismatch
        ):
            GM_PT.Printer.warning(
                "Warning! The current calculation makes use of both the "
                "AmideBB and AmideSC maps, but they try to emulate different "
                "versions. For most physical applications, this does "
                "not make sense. If you are absolutely sure that you want "
                "the two maps to follow different methods, make sure to set "
                "AmideBB.allow_map_mismatch to true. When in doubt, "
                "consult the README.",
                "map_AmideBB_2", True
            )

    # Initialize the prepro data structures (those that GMAP did for
    # non-prepro groups)
    MC_PC.initialize_prepro_properties(map_)

    # assign correct dipole function
    if map_.RunPars.dipole_map_choice == "Torii":
        map_.code.GM_calculate_dipole = MC_CM.calc_dipole_Torii
        map_.Core.dipole_gas_phase = np.float32(map_.Core.dipole_gas_phase)
        map_.Core.dipole_Torii_angle = np.float32(
            1 / np.tan(GM_con.deg2rad * map_.RunPars.Torii_dipole_angle))
    else:
        map_.code.GM_calculate_dipole = MC_CM.calc_dipole_Jansen

    if map_.RunPars.legacy_mode == "AIM":
        map_.code.GM_get_position_DMF = map_.code.GM_get_position
        map_.code.GM_get_position = MC_CM.get_position

    # now, knowing neighbors, we can determine the local atoms.
    oscillator_list = system.oscillators_ordered["AmideBB"]
    MC_LAF.find_local_atoms(map_, system, oscillator_list)

    MC_NM.read_maps(map_)
    MC_CM.determine_maps(oscillator_list, map_, system)

    if not map_.success:
        GM_PT.Printer.warning(
            "An issue occurred while initializing the AmideBB map stored at "
            f"{map_.directory}. Please first try restarting, then "
            "reinstalling, then contacting the map developer, as this map "
            "cannot be used like this. See the error above for more "
            "information. Quitting!",
            "map_AmideBB_0", True
        )


def GM_str_osc(map_, system, oscillator):
    """Explains how an oscillator should be printed.

    Example print: 'binding the residues GLY36 and LYS37'

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    oscillator : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The specific oscillator for which the string is required.

    Returns
    -------
    string : str
        The string that should be printed.
    """

    at0 = oscillator.used_atoms[0]
    at3 = oscillator.used_atoms[3]
    return (
        # example: binding the residues GLY36 and LYS37
        f"binding the residues {system.resnames[at0]}{system.resnums[at0]}"
        f" and {system.resnames[at3]}{system.resnums[at3]}"
    )


def GM_calculate_frequency(map_, system, osc):
    """Calculates the oscillating frequency for a given oscillator

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
        The specific oscillator for which the frequency is required.

    Returns
    -------
    freq : float
        The frequency calculated for this oscillator this frame.
    """

    if osc.resnames[1] == "PRO":
        gasfreq = map_.Core.frequency_gas_phase_prepro
        freqarr = map_.Core.frequency_data_array_linear_prepro
    else:
        gasfreq = map_.Core.frequency_gas_phase
        freqarr = map_.Core.frequency_data_array_linear

    # np.seterr(all='raise')
    # try:
    freq = gasfreq + np.sum(np.multiply(osc.VEGout, freqarr))
    # except Exception as ex:
    #     dpr(osc.oscix)
    #     dpr(osc.VEGout)
    #     dpr(freqarr)
    #     if osc.oscix > 10:
    #         raise ex
    #     else:
    #         freq = gasfreq

    if (
        map_.RunPars.frequency_map_choice != "Tokmakoff"
        and map_.RunPars.consider_nearest_neighbours
    ):
        freq += MC_CM.neighbor_influence(map_, system, osc)

    return freq


def GM_calculate_raman(map_, system, osc):
    """Returns the raman tensor as a length-6 vector: (xx, xy, xz, yy, yz, zz)

    This method can easily be adapted by other maps for working with raman
    tensors. Make sure that rotation_matrix is an orthogonal 3*3 numpy array
    (so, the vectors making it up are orthonormal).
    Then, the raman_tensor_local can be freely chosen.

    As this is so easily adaptable, this could also be made a standard
    built-in function for GMAP. The only reason this is not the case
    currently, is because raman tensors from maps like this are not
    common yet, so a 'usual' way of determining them has not yet been
    created. Maybe, they won't stay of fixed magnitude in local coordinates
    forever, but depend on sth like VEG or atomic distances in the future.

    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    osc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The specific oscillator for which the transformation is required.

    Returns
    -------
    raman_tensor : `np.ndarray`
        The raman tensor calculated for this oscillator this frame.
    """

    # the rotation matrix is available as long as the map specifies
    # estatic_choice to be E or G (which is the case here). It is made
    # available immediately at the beginning of the frame.
    COvec = osc.rotation_matrix[0, :]
    CNvec = osc.rotation_matrix[1, :]
    Zvec = osc.rotation_matrix[2, :]

    theta = 34*np.pi/180
    raman_tensor_local = np.diag([20, 4, 1])

    rotation_matrix = np.zeros((3, 3))
    rotation_matrix[0] = np.cos(theta) * COvec - np.sin(theta) * CNvec
    rotation_matrix[1] = np.sin(theta) * COvec + np.cos(theta) * CNvec
    rotation_matrix[2] = Zvec

    raman_tensor_system = (
        rotation_matrix.T @ raman_tensor_local @ rotation_matrix)

    # old (AIM) version:
    # def tp(vect1):  # tensor product
    #     tensor = np.zeros((6), dtype='float32')
    #     tensor[:3] = vect1[0]*vect1
    #     tensor[3:5] = vect1[1]*vect1[1:]
    #     tensor[5] = vect1[2]*vect1[2]
    #     return tensor
    # Rvec = tp(Rtens[0]) * 20 + tp(Rtens[1]) * 4  + tp(Rtens[2])
    # (here, Rtens is what the current version calls rotation_matrix)

    # now, to numpify this, first, redefine tp.
    # def tp(vect1):
    #     return (vect1[:, None] * vect1[None, :])[np.triu_indices(3)]

    # then, we can do the entire array at once:
    # consts = np.array([20, 4, 1])
    # Rvec = (
    #     Rtens[:, :, None] * Rtens[:, None, :] * consts[:, None, None]
    # ).sum(axis=0)[np.triu_indices(3)]

    # in summation notation (forgetting the triu-indices for flattening):
    # with A_ij == A[i, j]
    # Rvec[i, j] = sum{k=1 -> k=3}(Rtens[k, i] * Rtens[k, j] * consts[k])

    # now, is this equivalent to the new method? Lets derive the summation
    # notation for the new method! (R = rotation matrix, A = local raman tens)
    # Assuming A is diagonal (so only a[i, i] exist)
    # Rvec = R.T @ A @ R
    # Rvec[i, j] = sum{k=1 -> k=3}(R.T[i, k] * (A @ R)[k, j])
    #            = sum{k=1 -> k=3}(R[k, i] * A[k, k] * R[k, j])
    # this is the same as the summation for the AIM version!

    # footnote: what is (A @ R)[k, j]?
    # write it out: (A @ R)[i, j] = sum{k=1 -> k=3}(A[i, k] * R[k, j])
    # but, as only k==i exists for A (rest is 0), this becomes:
    # (A @ R)[i, j] = A[i, i] * R[i, j]

    return raman_tensor_system[np.triu_indices(3)]
