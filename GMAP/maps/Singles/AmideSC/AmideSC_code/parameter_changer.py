
def adjust_map_core_raw(Files, map_):
    # atom order: CG  OD1  CB  ND2  HD21  HD22

    # Based on frequency map choice (Tokmakoff, Skinner, Jansen, Cho, Hirst),
    # some keywords from rawpars need to be set to different values.
    adjust_mcr_freqchoice(map_)

    # Based on position choice, set the position keywords from the map to
    # different values
    choice = map_.RunPars.pos_choice
    match choice:
        case "C":  # the default
            map_.rawcore["position"] = ["0"]
        case "O":
            map_.rawcore["position"] = ["1"]
        case "N":
            map_.rawcore["position"] = ["3"]
        case "D":
            map_.rawcore["position"] = ["4"]
        case "CNO":
            # map_.rawcore["position"] = ["(0+0.665*1+0.258*3)/1.923"]
            map_.rawcore["position"] = ["0.077*0+0.665*1+0.258*3"]

    match map_.RunPars.dipole_map_choice:
        case "Torii":
            map_.rawcore["dipole_data_file"] = ["[N/A]"]
            map_.rawcore["dipole_gas_phase"] = ["0.276"]
        case "Jansen":
            if not map_.RunPars.frequency_map_choice == "Jansen":
                # throw up an error, and set map_.success to false.
                pass
    # else: default (xyz jansen_dipoles.txt)


def adjust_mcr_freqchoice(map_):
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
            map_.rawcore["frequency_gas_phase"] = ["1740"]
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
