
# 3rd party lib imports
import numpy as np

# gmap imports
from GMAP.src.tools import constants as GM_con
from GMAP.src.tools import FileHandler as GM_FH
from GMAP.src.tools import ParameterParser as GM_PP
from GMAP.src.tools import PrintTools as GM_PT

# own module imports
import TRESP_code.TRESPclib as MC_TC

_ = GM_con.bohr  # to validify the import. The import is needed for exec.


def GM_change_coup_type(map_, system, oscix1, osc1, oscix2, osc2):
    return map_.name  # return TRESP


def GM_post_init(map_, system):
    # load in the charges from all different kinds of oscillators
    map_.charges = {}
    map_charges = {}
    for osc in system.oscillators_ordered_coup[map_.name]:
        # if an issue occurs, no need to continue
        if not map_.success:
            return

        # If a map has a dedicated function, use that instead of interpreting
        # the provided file.

        if hasattr(osc.Map.code, "CP_TRESP_get_charges"):
            map_.charges[osc.oscix] = osc.Map.code.CP_TRESP_get_charges(
                    osc.Map, system, osc)
        # only look for each type of singles once.
        elif osc.Map.name in map_charges:
            map_.charges[osc.oscix] = map_charges[osc.Map.name]
        else:
            map_charges[osc.Map.name] = get_charges(map_, osc.Map)
            print(map_charges[osc.Map.name])
            map_.charges[osc.oscix] = map_charges[osc.Map.name]

    MC_TC.init_map_for_clib(map_, system)


def GM_pre_run(map_, system):
    map_.osclens = [
        len(map_.charges.get(osc.oscix, [])) for osc in system.oscillators]
    # do this before converting osclens to array (faster)
    map_.all_used_ats = []
    for osclen, osc in zip(map_.osclens, system.oscillators):
        map_.all_used_ats.extend(osc.used_atoms[:osclen])

    map_.osclens = np.array(map_.osclens, dtype="int32")
    map_.osclens_c = np.ctypeslib.as_ctypes(map_.osclens)
    map_.oscstart = np.concatenate(
        (np.zeros(1, dtype="int32"), np.cumsum(map_.osclens)[:-1]),
        dtype="int32")  # set the dtype again, as linux changes it here.
    map_.oscstart_c = np.ctypeslib.as_ctypes(map_.oscstart)
    map_.all_used_ats = np.array(map_.all_used_ats, dtype="int32")
    map_.all_used_ats_c = np.ctypeslib.as_ctypes(map_.all_used_ats)

    map_.charge_array = []
    for osc in system.oscillators:
        map_.charge_array.extend(map_.charges.get(osc.oscix, []))
    map_.charge_array = np.array(map_.charge_array, dtype="float32")
    map_.charge_array_c = np.ctypeslib.as_ctypes(map_.charge_array)

    map_.allpairs = np.array(map_.allpairs, dtype="int32")
    map_.allpairs_c = np.ctypeslib.as_ctypes(np.ravel(map_.allpairs))
    map_.n_allpairs = np.int32(map_.allpairs.shape[0])


def GM_calc_coupling(map_, system, hamiltonian):
    hamiltonian_c = np.ctypeslib.as_ctypes(np.ravel(hamiltonian))
    map_.clib.calc_coupling(map_, system, hamiltonian_c)


def get_charges(map_, oscmap):
    fnameraw = oscmap.rawcore[f"{map_.name}.charges_filename"][0]
    fname = (oscmap.directory / fnameraw).resolve()
    if GM_FH.try_file(fname) is None:
        GM_PT.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file doesn't exist:\n"
            f"{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TRESP_1", False
        )
        map_.success = False
        return None

    if not GM_FH.check_file_readability(fname, False, False):
        GM_PT.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file is of the wrong format:"
            f"\n{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TRESP_2", False
        )
        map_.success = False
        return None

    contents = []
    with open(fname, encoding="utf-8") as fhand:
        for line in fhand:
            line = GM_PP.cleanline(line).strip()
            if line:
                contents.append(line)

    try:
        contents = [float(item) for item in contents]
    except Exception:
        GM_PT.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file has the wrong contents:"
            f"\n{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TRESP_3", False
        )
        map_.success = False
        return None

    if len(contents) > len(oscmap.Core.used_atoms):
        GM_PT.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file has too many contents:"
            f"\n{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TRESP_4", False
        )
        map_.success = False
        return None

    contents = np.array(contents)

    keyword = f"{map_.name}.charges_multiply"
    # optional multiplication - if not needed, skip.
    if keyword not in oscmap.rawcore:
        return contents

    # ----------------------------------------------------

    cmdstr = "multiplier = " + " ".join(oscmap.rawcore[keyword])
    mapdir = oscmap.directory
    pars = {}
    try:
        exec(cmdstr, globals(), pars)
    except Exception:
        GM_PT.Printer.warning(
            "\nCould not interpret the choice for the keyword "
            f"'{keyword}' in the file {mapdir / 'core.txt'}. "
            "Please make sure the choice only contains numbers (and "
            "optionally a single '.') that represent a decimal value. "
            "Alternatively, make sure it is a python-parsable string. ",
            "map_TRESP_5", False
        )
        map_.success = False
        return None

    try:
        multiplier = float(pars["multiplier"])
    except Exception:
        GM_PT.Printer.warning(
            "\nCould not interpret the choice for the keyword "
            f"'multiply_freq' in the file {mapdir / 'core.txt'}. "
            "Please make sure the choice only contains numbers (and "
            "optionally a single '.') that represent a decimal value. "
            "Alternatively, make sure it is a python-parsable string.2 ",
            "map_TRESP_5", False
        )
        map_.success = False
        return None

    return contents * multiplier
