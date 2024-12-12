
# standard library imports
from pathlib import Path

# 3rd party imports
import numpy as np

# GMAP imports
import GMAP.src.tools.constants as GM_Con
import GMAP.src.tools.PrintTools as GM_PT


def adjust_map_core_raw(map_):
    """Changes some choices in core.txt to match the chosen maps

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    # atom order: CG  OD1  CB  ND2  HD21  HD22

    # Based on frequency map choice (Tokmakoff, Skinner, Jansen, Cho, Hirst),
    # some keywords from rawpars need to be set to different values.
    adjust_mcr_freqchoice(map_)

    # Based on position choice, set the position keywords from the map to
    # different values
    pardict = {
        "pos_choice": "position",
        "dp1_choice": "doublepos_0",
        "dp2_choice": "doublepos_1"
    }
    for parameter in ("pos_choice", "dp1_choice", "dp2_choice"):
        choice = getattr(map_.RunPars, parameter)
        rawpar = pardict[parameter]
        match choice:
            case "C":  # the default
                map_.rawcore[rawpar] = ["0"]
            case "O":
                map_.rawcore[rawpar] = ["1"]
            case "N":
                map_.rawcore[rawpar] = ["3"]
            case "D":
                map_.rawcore[rawpar] = ["4"]
            case "Torii":
                # map_.rawcore[rawpar] = ["(0+0.665*1+0.258*3)/1.923"]
                map_.rawcore[rawpar] = ["0.077*0+0.665*1+0.258*3"]

    match map_.RunPars.dipole_map_choice:
        case "Torii":
            map_.rawcore["dipole_data_file"] = ["[N/A]"]
            map_.rawcore["dipole_gas_phase"] = ["0.276"]
        case "Jansen":
            if not map_.RunPars.frequency_map_choice == "Jansen":
                GM_PT.Printer.warning(  # no exitbool - error is not fatal.
                    "Error in the map AmideSC: The Jansen dipole map was "
                    "requested without using the Jansen frequency map. Either "
                    "change your frequency map choice to Jansen, or "
                    "use a different dipole map.",
                    "map_AmideBB_1"
                )
                # Indicate there is an error with this map, so GMAP shouldn't
                # use it (GMAP raises fatal error if this map is requested).
                map_.success = False
    # else: default (xyz jansen_dipoles.txt)


def adjust_mcr_freqchoice(map_):
    """Makes changes to core.txt based on frequency map choice.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    # atom order: CG  OD1  CB  ND2  HD21  HD22

    choice = map_.RunPars.frequency_map_choice
    parname = "frequency_data_file_linear"
    secpar = "frequency_data_file_linear_prepro"
    match choice:
        case "Skinner":  # the default
            map_.rawcore["electrostatic_atoms"] = ["0", "3"]  # C and N
            map_.rawcore["electrostatic_choice"] = ["E"]
            map_.rawcore["frequency_gas_phase"] = ["1684"]
            map_.rawcore["frequency_gas_phase_prepro"] = ["1657"]
            map_.rawcore[parname] = ["frequency_maps/Skinner.txt"]
            map_.rawcore[secpar] = ["frequency_maps/Skinner.txt"]
        case "Tokmakoff":
            map_.rawcore["electrostatic_atoms"] = ["1"]  # Oxygen!
            map_.rawcore["electrostatic_choice"] = ["E"]
            map_.rawcore["frequency_gas_phase"] = ["1710.3"]
            map_.rawcore["frequency_gas_phase_prepro"] = ["1666"]
            map_.rawcore[parname] = ["frequency_maps/Tokmakoff_gen.txt"]
            map_.rawcore[secpar] = ["frequency_maps/Tokmakoff_prepro.txt"]
        case "Jansen":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["G"]
            map_.rawcore["frequency_gas_phase"] = ["1717"]
            map_.rawcore["frequency_gas_phase_prepro"] = ["1690"]
            map_.rawcore[parname] = ["frequency_maps/Jansen_gen.txt"]
            map_.rawcore[secpar] = ["frequency_maps/Jansen_prepro.txt"]
        case "Cho":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["V"]
            map_.rawcore["frequency_gas_phase"] = ["1717"]
            map_.rawcore["frequency_gas_phase_prepro"] = ["1690"]
            map_.rawcore[parname] = ["frequency_maps/Cho.txt"]
            map_.rawcore[secpar] = ["frequency_maps/Cho.txt"]
        case "Hirst":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["V"]
            map_.rawcore["frequency_gas_phase"] = ["1717"]
            map_.rawcore["frequency_gas_phase_prepro"] = ["1690"]
            map_.rawcore[parname] = ["frequency_maps/Hirst.txt"]
            map_.rawcore[secpar] = ["frequency_maps/Hirst.txt"]


def oscillator_sorter(map_, system, oscillator_list):
    """Sorts the oscillators into the desired order.

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

    Returns
    -------
    newlist : list of :class:`~GMAP.src.tools.SystemReader.Oscillator`
        The same oscillators as in oscillator list, but in the desired
        order.
    """

    # first, label oscillators
    # used ats order:    res0{C O CA} res1{N H CA}
    for oscillator in oscillator_list:
        oscillator.resnames = (
            system.resnames[oscillator.used_atoms[0]],
            system.resnames[oscillator.used_atoms[3]]
        )
        oscillator.resnums = (
            system.resnums[oscillator.used_atoms[0]],
            system.resnums[oscillator.used_atoms[3]]
        )

    # then, sort oscillators
    newlist = []
    if map_.RunPars.residue_order == "resname":
        all_amino_acid_codes = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LYS", "LEU", "MET",
            "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        ]
        for aa1 in all_amino_acid_codes:
            for aa2 in all_amino_acid_codes:
                for ix, oscillator in enumerate(oscillator_list):
                    if oscillator.resnames == (aa1, aa2):
                        newlist.append(oscillator_list.pop(ix))
    else:
        # now, choice is 'resnum'. To change order to AIM order:
        while len(oscillator_list) > 0:
            smallest_ix = 0
            smallest_resix = 999999999
            for ix, oscillator in enumerate(oscillator_list):
                resix = system.resnums[oscillator.used_atoms[0]]
                if resix < smallest_resix:
                    smallest_resix = resix
                    smallest_ix = ix
            newlist.append(oscillator_list.pop(smallest_ix))

    return newlist


def initialize_prepro_properties(map_):
    """Read and parse the data from prepro map files.

    Normally, GMAP automatically parses the data files, but this map
    needs multiple, so these extras need to be done manually.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    # For assigning constants, first they must be created
    try:
        pp_gasfreq = np.float32(map_.rawcore["frequency_gas_phase_prepro"][0])
    except Exception as ex:  # no 0th entry, not floatable
        mapdir = map_.directory
        GM_PT.Printer.warning(
            "\nCould not interpret the choice for the parameter "
            "'frequency_gas_phase_prepro'"
            f" in the file {mapdir / 'core.txt'}. Please make sure "
            "the choice consists of a single decimal number.",
            "MI_MC_7", exception=ex
        )
        map_.success = False
    else:
        map_.Core.frequency_gas_phase_prepro = pp_gasfreq

    pp_freqarr = map_.Core.parse_frequency_data_file(
        map_.rawcore, map_.directory, "frequency_data_file_linear_prepro")
    if pp_freqarr is None:
        mapdir = map_.directory
        GM_PT.Printer.warning(
            "\nCould not interpret the choice for the parameter "
            "'frequency_data_file_linear_prepro'"
            f" in the file {mapdir / 'core.txt'}. Please make sure "
            "the choice consists of a single decimal number.",
            "MI_MC_7"
        )
        map_.success = False
    else:
        if map_.Core.length_units == "bohr":
            conv_factor = GM_Con.bohr2ang
            pp_freqarr[:, 0] *= conv_factor
            pp_freqarr[:, 1:4] *= conv_factor**2
            pp_freqarr[:, 4:] *= conv_factor**3
        map_.Core.frequency_data_array_linear_prepro = pp_freqarr

    map_.Core.dipole_gas_phase_prepro = [np.float32(item) for item in [
        -0.268549, 0.086947, 0.0]]
    map_.Core.dipole_gas_phase_array_prepro = np.array(
        map_.Core.dipole_gas_phase_prepro)
    try:
        fname = "jansen_dipoles_prepro.txt"
        fname = (Path(__file__).resolve().parent.parent / fname).resolve()
        fdata = np.genfromtxt(fname, "float32", missing_values=0, ndmin=2)
    except Exception as ex:
        GM_PT.Printer.warning(
            "\nNumpy could not interpret the contents of the file "
            f"{fname}. Please make sure the file contains only decimal "
            "numbers in a grid.",
            "MI_MC_7", exception=ex
        )
        map_.success = False
        return

    dip_arr = map_.Core.confirm_array_size(fdata, 10, 12, fname)
    if not map_.success:
        return

    map_.Core.dipole_data_array_prepro = dip_arr.reshape((3, -1, 10))

    if map_.Core.length_units == "bohr":
        conv_factor = GM_Con.bohr2ang
        map_.Core.dipole_data_array_prepro[:, :, 0] *= conv_factor
        map_.Core.dipole_data_array_prepro[:, :, 1:4] *= conv_factor**2
        map_.Core.dipole_data_array_prepro[:, :, 4:] *= conv_factor**3
