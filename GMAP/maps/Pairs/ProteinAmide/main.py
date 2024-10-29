
def GM_change_coup_type(map_, system, oscix1, osc1, oscix2, osc2):
    mapnames = {
        "None": None,
        "DipDip": "DipDip",
        "TDCKrimm": "ProteinAmide_TDCKrimm",
        "TDCTasumi": "ProteinAmide_TDCTasumi",
        "TCC": "ProteinAmide_TCC",
        "Tasumi": "ProteinAmide_Tasumi",
        "GLDP": "ProteinAmide_GLDP"
    }

    # If any of the two is a sidechain, they cannot be nearest neighbours
    if osc1.Map.name == "AmideSC" or osc2.Map.name == "AmideSC":
        return mapnames[map_.RunPars.coupling_choice]

    # Both are backbone, but not neighbors
    if osc1 not in (osc2.NtermNB, osc2.CtermNB):
        return mapnames[map_.RunPars.coupling_choice]
    else:
        return mapnames[map_.RunPars.NN_coupling_choice]
