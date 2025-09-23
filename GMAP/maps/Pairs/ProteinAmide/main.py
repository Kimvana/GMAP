
def GM_change_coup_type(map_, system, oscix1, osc1, oscix2, osc2):
    mapnames = {
        "None": None,
        "DipDip": "DipDip",
        "TDCKnoester": "ProteinAmide_TDCKnoester",
        "TDCTasumi": "ProteinAmide_TDCTasumi",
        "TCC": "ProteinAmide_TCC",
        "Tasumi": "ProteinAmide_Tasumi",
        "GLDP": "ProteinAmide_GLDP"
    }

    legacy_mode = map_.run_pars.legacy_mode

    # If any of the two is a sidechain, they cannot be nearest neighbours
    if osc1.map.name == "AmideSC" or osc2.map.name == "AmideSC":
        # original code forced the double-sidechain couplings to be
        # coupled by a TDC method (TDCTasumi by default, but if TDCKnoester is
        # already active, use that)
        # IMPORTANT! The old version did not have DipDip or none, so this
        # if check is not designed to make those run smoothly.
        if (
            legacy_mode == "AmideImaps"
            and osc1.map.name == "AmideSC"
            and osc2.map.name == "AmideSC"
            and map_.run_pars.coupling_choice == "TCC"
        ):
            return "ProteinAmide_TDCTasumi"
        return mapnames[map_.run_pars.coupling_choice]

    # Both are backbone, but not neighbors
    if osc1 not in (osc2.NtermNB, osc2.CtermNB):
        return mapnames[map_.run_pars.coupling_choice]

    # From here, they are neighbors

    # AmideImaps had a bug when two opposite ends of a cyclic chain connect:
    # Then, always use GLDP (used to be the NN default at time of writing)
    elif (
        legacy_mode == "AmideImaps"
        and osc1.resnums[0] not in (osc2.resnums[0] - 1, osc2.resnums[0] + 1)
    ):
        return "ProteinAmide_GLDP"
    else:
        return mapnames[map_.run_pars.NN_coupling_choice]
