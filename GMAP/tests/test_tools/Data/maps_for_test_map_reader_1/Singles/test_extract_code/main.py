

def GM_adjust_RunPars(map_):
    """Makes the necessary changes to Map.RunPar.

    Is expected to not return anything - return value is not caught.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    if map_.RunPars.overwrite_neutral_charge_threshold != 0:
        setattr(
            map_.RunPars.MainRunPars,
            "neutral_charge_threshold",
            map_.RunPars.overwrite_neutral_charge_threshold
        )


def for_testing(number):
    return 10 - number
