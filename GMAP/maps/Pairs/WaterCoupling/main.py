# Not sure about the options below, intramolecular must always be done with the
# intramolecular map and inter with dipdip
def GM_change_coup_type(map_, system, oscix1, osc1, oscix2, osc2):
    mapnames = {
        "None": None,
        "DipDip": "DipDip",
        "WaterCoupling": "WaterCoupling",
        "Water_Intra": "Water_Intra"
    }

    # If the oscillators are in the same water molecule we must choose intra
    # unless we do not want coupling
    # Check if the OH stretches share the same oxygen
    if osc1.used_atoms[0] == osc2.used_atoms[0]:
        return mapnames[map_.run_pars.intra_coupling_choice]

    # The OH stretches are in different molecules. If the user specified to use
    # intramolecular couplings we force them to use the dipole dipole instead
    else:
        return mapnames[map_.run_pars.inter_coupling_choice]
