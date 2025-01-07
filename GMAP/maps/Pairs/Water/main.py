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
        return mapnames[map_.RunPars.intra_coupling_choice]

    # The OH stretches are in different molecules. If the user specified to use
    # intramolecular couplings we force them to use the dipole dipole instead
    else: 
        return mapnames[map_.RunPars.inter_coupling_choice]


