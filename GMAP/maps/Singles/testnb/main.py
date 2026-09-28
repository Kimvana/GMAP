from GMAP.src.tools.print_tools import devprint as dpr

MAX_MODES = 9


def GM_adjust_oscillators(map_, system, oscillator_list):
    """Makes the necessary changes to the list of oscillators.

    The program finds all oscillators mathing the instructions from
    core.txt. However, there is no way for the program to avoid double
    counting symmetrical groups (like the cystbridge mockup example).
    If a map knows its group is symmetrical, this function can be
    designed to only return half of the inputs.

    Another possible use is for the code of the map to get to know its
    oscillators. When all oscillators are passed through this function,
    the (global) atom number of the first atom of this group (for
    example) can be linked to a specific property the group might need
    to know. This might be useful if a map needs to cover two very
    similar oscillators.

    .. note::
        This function is called separately for each struct that the map
        defines. So take into account that the function could be called
        multiple times within a single simulation!

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    oscillator_list : list of :class:`~GMAP.src.tools.system_reader.Oscillator`
        All oscillators belonging to a single struct of this map.

    Returns
    -------
    oscillator_list : list of :class:`~GMAP.src.tools.system_reader.Oscillator`
        All oscillators belonging to a single struct of this map.
    """

    for oscillator in oscillator_list:
        dpr(oscillator)
        dpr(oscillator.used_atoms)

    return oscillator_list


def GM_filter_oscillators(map_, system, oscillators):
    """Applies requested black/whitelist filters.

    This function is the same as the default implementation, but uses a
    different implementation of the 'filter_single_line' function.

    Supports the default options provided by GMAP, as well as these:

    """

    runpars = map_.run_pars.main_run_pars
    wl_rules = runpars.singles_whitelist_dict.get(map_.name, [[":All"]])
    bl_rules = runpars.singles_blacklist_dict.get(map_.name, [[":None"]])

    # first, select all in whitelist from oscillators
    filtered = set()
    oscset = set(oscillators)
    for rule in wl_rules:
        success, filtered = filter_single_line(
            rule, "white", filtered, oscset, map_, system)
        if not success:
            GM_pt.Printer.warning(
                "\nUsing the parameter 'singles_whitelist', the map "
                f"{map_.name} "
                "requestested a specific selection, but it was not "
                "recognized. ",
                "map_lys-glu-dimer-eq_0", True,
                GMAPerrclass=GM_ex.GmapFileSyntaxError
            )

    # from the whitelisted, remove all those in blacklist.
    for rule in bl_rules:
        success, filtered = filter_single_line(
            rule, "black", filtered, oscset, map_, system
        )
        if not success:
            GM_pt.Printer.warning(
                "\nUsing the parameter 'singles_whitelist', the map "
                f"{map_.name} "
                "requestested a specific selection, but it was not "
                "recognized. ",
                "map_lys-glu-dimer-eq_0", True,
                GMAPerrclass=GM_ex.GmapFileSyntaxError
            )

    # sort the oscillators in correct order (same as before)
    filtered_list = [osc for osc in oscillators if osc in filtered]
    return filtered_list


def filter_single_line(line, BW, found, avail, map_, system):
    """Applies a single requested black/whitelist filter.

    This function is very similar to the default implementation, but
    treats the resnums differently, and does not support resnames.

    Supports the default options provided by GMAP, as well as these:

    resnums
        Same name as default GMAP, but different implementation.
        Expected behaviour:

        >>> resnums 1:5 3:2
        selects 2 dimers: LYS1-GLU5 and LYS3-GLU2

        >>> resnums 1-3:5
        selects 3 dimers: LYS1-GLU5, LYS2-GLU5 and LYS3-GLU5

        >>> renums 1,3:5,2 8:9
        selects 5 dimers: LYS1-GLU5, LYS1-GLU2, LYS3-GLU5, LYS3-GLU2,
        and LYS8-GLU9

        >>> resnums 1-3,7:5,6
        select 8 dimers: LYS1-GLU5, LYS2-GLU5, LYS3-GLU5, LYS7-GLU5,
        LYS1-GL6, LYS2-GLU6, LYS3-GLU6, LYS7-GLU6

    excitations
        Same default GMAP resnum implementation

    singles_whitelist lys-glu-dimer-eq [keyword] [choice]
    """

    match line[0].lower():
        case ":all":
            if BW == "white":
                return True, avail.copy()
            else:
                return True, set()
        case ":none":
            if BW == "white":
                return True, set()
            else:
                return True, found.copy()
        case "resnums":
            # Check whether the choice is alphanumeric
            # we can use/support hyphens, too, but not commas/periods.
            if not set("".join(line[1:])).issubset("1234567890-,:"):
                GM_pt.Printer.warning(
                    f"\nUsing the parameter 'singles_{BW}list', the map "
                    f"{map_.name}"
                    "was requestested certain residue numbers, but this "
                    "specification used non-numeric characters. Please make "
                    "sure to only use numbers and hyphens. ",
                    "map_lys-glu-dimer-eq_0", True,
                    GMAPerrclass=GM_ex.GmapFileSyntaxError
                )

            # find all allowed pairs
            all_requested = set()
            for selection in line[1:]:
                terms = selection.split(":")
                if len(terms) != 2:
                    GM_pt.Printer.warning(
                        f"\nUsing the parameter 'singles_{BW}list', the map "
                        f"{map_.name}"
                        "was requestested certain residue numbers, but this "
                        "specification is of invalid format. Make sure to use "
                        "exactly one colon per selection group.",
                        "map_lys-glu-dimer-eq_0", True,
                        GMAPerrclass=GM_ex.GmapFileSyntaxError
                    )
                try:
                    parsed = [
                        map_.core.allow_ranges(term.split(","), system.nres)
                        for term in terms]
                except IndexError as ierr:
                    GM_pt.Printer.warning(
                        f"\nUsing the parameter 'singles_{BW}list', the map "
                        f"{map_.name}"
                        "was requestested certain residue numbers, but the "
                        "specific residue numbers requested do not exist in "
                        " the provided MD system. ",
                        "map_lys-glu-dimer-eq_0", True, exception=ierr,
                        GMAPerrclass=GM_ex.GmapIndexError
                    )
                except Exception as ex:
                    GM_pt.Printer.warning(
                        f"\nUsing the parameter 'singles_{BW}list', the map "
                        f"{map_.name}"
                        "was requestested certain residue numbers, but the "
                        "specific choice provided could not be interpreted. "
                        "Please make sure the choice consists of nothing but "
                        "numbers separated by spaces "
                        "and/or ranges of integers separated by a hyphen.",
                        "map_lys-glu-dimer-eq_0", True, exception=ex,
                        GMAPerrclass=GM_ex.GmapFileSyntaxError
                    )

                # parsed is a list of length 2. item 1 is all possible indices
                # for lysine, item 2 all for glutamate.
                for ix1 in parsed[0]:
                    for ix2 in parsed[1]:
                        all_requested.add((ix1, ix2))

            # only keep an oscillator if its number is expected.
            filtered = {
                osc for osc in avail
                if (
                    system.resnums[osc.used_atoms[6]],
                    system.resnums[osc.used_atoms[16]]
                ) in all_requested
            }
            if BW == "white":
                found |= filtered
            else:
                found -= filtered
            return True, found

        case "excitations":
            # Check whether the choice is alphanumeric
            # we can use/support hyphens, too, but not commas/periods.
            if not set("".join(line[1:])).issubset("1234567890-"):
                GM_pt.Printer.warning(
                    f"\nUsing the parameter 'singles_{BW}list', the map "
                    f"{map_.name}"
                    "was requestested certain excitations, but this "
                    "specification used non-numeric characters. Please make "
                    "sure to only use numbers and hyphens. ",
                    "map_lys-glu-dimer-eq_0", True,
                    GMAPerrclass=GM_ex.GmapFileSyntaxError
                )
            # find all allowed numbers
            try:
                modes = set(map_.core.allow_ranges(line[1:], MAX_MODES))
            except IndexError as ierr:
                GM_pt.Printer.warning(
                    f"\nUsing the parameter 'singles_{BW}list', the map "
                    f"{map_.name}"
                    "was requestested certain excitations, but the "
                    "specific excitations requested do not exist. ",
                    "map_lys-glu-dimer-eq_1", True, exception=ierr,
                    GMAPerrclass=GM_ex.GmapIndexError
                )
            except Exception as ex:
                GM_pt.Printer.warning(
                    f"\nUsing the parameter 'singles_{BW}list', the map "
                    f"{map_.name}"
                    "was requestested certain excitations, but the "
                    "specific choice provided could not be interpreted. "
                    "Please make sure the choice consists of nothing but "
                    "numbers separated by spaces "
                    "and/or ranges of integers separated by a hyphen.",
                    "map_lys-glu-dimer-eq_2", True, exception=ex,
                    GMAPerrclass=GM_ex.GmapFileSyntaxError
                )

            # only keep an oscillator if its number is expected.
            filtered = {
                osc for osc in avail
                if osc.mode in modes
            }
            if BW == "white":
                found |= filtered
            else:
                found -= filtered
            return True, found

        case _:
            return False, found
