# Not sure about the options below, intramolecular must always be done with the intramolecular map and inter with dipdip
def GM_change_coup_type(map_, system, oscix1, osc1, oscix2, osc2):
    mapnames = {
        "None": None,
        "DipDip": "DipDip",
        "Water": "WaterIntra",
    }

    legacy_mode = map_.RunPars.legacy_mode

    # If the oscillators are in the same water molecule we must choose intra
    # unless we do not want coupling
    if osc1.used_atoms[1] == osc2.used_atoms[1]:
        if map_.RunPars.coupling_choice != "None":
            return "WaterIntra"
        return mapnames[map_.RunPars.coupling_choice]

    # The OH stretches are in different molecules. If the user specified to use
    # intramolecular couplings we force them to use the dipole dipole instead
    if map_.RunPars.coupling_choice == "WaterIntra":
            return "DipDip"
    return mapnames[map_.RunPars.coupling_choice]

    # Use the choice specified by the user (should be DipDip here, but more
    # advanced maps could be added in the future.)
    return mapnames[map_.RunPars.coupling_choice]

