

def GM_adjust_run_pars(map_):
    """Makes the necessary changes to map_.run_pars.

    Is expected to not return anything - return value is not caught.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    if map_.run_pars.overwrite_neutral_charge_threshold != 0:
        setattr(
            map_.run_pars.main_run_pars,
            "neutral_charge_threshold",
            map_.run_pars.overwrite_neutral_charge_threshold
        )


def for_testing(number):
    return 10 - number
