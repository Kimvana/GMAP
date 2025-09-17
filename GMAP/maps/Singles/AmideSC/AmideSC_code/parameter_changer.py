
# GMAP imports
import GMAP.src.tools.print_tools as GM_pt


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
                GM_pt.Printer.warning(  # no exitbool - error is not fatal.
                    "Error in the map AmideSC: The Jansen dipole map was "
                    "requested without using the Jansen frequency map. Either "
                    "change your frequency map choice to Jansen, or "
                    "use a different dipole map.",
                    "map_AmideSC_1"
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
    match choice:
        case "Skinner":  # the default
            map_.rawcore["electrostatic_atoms"] = ["0", "3"]  # C and N
            map_.rawcore["electrostatic_choice"] = ["E"]
            map_.rawcore["frequency_gas_phase"] = ["1714"]
            map_.rawcore[parname] = ["frequency_maps/Skinner.txt"]
        case "Tokmakoff":
            map_.rawcore["electrostatic_atoms"] = ["1"]  # Oxygen!
            map_.rawcore["electrostatic_choice"] = ["E"]
            map_.rawcore["frequency_gas_phase"] = ["1740.3"]
            map_.rawcore[parname] = ["frequency_maps/Tokmakoff.txt"]
        case "Jansen":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["G"]
            map_.rawcore["frequency_gas_phase"] = ["1747"]
            map_.rawcore[parname] = ["frequency_maps/Jansen.txt"]
        case "Cho":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["V"]
            map_.rawcore["frequency_gas_phase"] = ["1747"]
            map_.rawcore[parname] = ["frequency_maps/Cho.txt"]
        case "Hirst":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["V"]
            map_.rawcore["frequency_gas_phase"] = ["1747"]
            map_.rawcore[parname] = ["frequency_maps/Hirst.txt"]
        case "Reppert_4PN-4":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["V"]
            map_.rawcore["frequency_gas_phase"] = ["1746.6"]
            map_.rawcore[parname] = ["frequency_maps/Reppert_4PN-4.txt"]
        case "Reppert_4PN-150":
            map_.rawcore["electrostatic_atoms"] = ["0", "1", "3", "4"]
            map_.rawcore["electrostatic_choice"] = ["V"]
            map_.rawcore["frequency_gas_phase"] = ["1776.4"]
            map_.rawcore[parname] = ["frequency_maps/Reppert_4PN-150.txt"]
