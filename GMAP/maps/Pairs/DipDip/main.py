"""Main code file for dipole-dipole map.

Also is the template/explanation for all that a map can hand directly to
GMAP. Of course, more functions are allowed, and these files can import
other custom python files stored in the same directory (or a directory
therein) as this file.
"""

# 3rd party lib imports
from numba import njit
import numpy as np

# gmap imports
from GMAP.src.tools import MathFunctions as GM_MF
from GMAP.src.tools import ReferenceHandler as GM_RH


# A function to adjust the parameters of the map. For some kinds of
# parameter (especially if theres multiple that are linked), the way
# RunPar is built might not be correct. In this function, the user can
# fix that.
def GM_adjust_RunPars(map_):
    """Makes the necessary changes to map_.RunPar.

    Is expected to not return anything - return value is not caught.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    pass


# A function to adjust the choices made in core.txt. Perhaps, based on
# a detected parameter, a different choice is preferred. This function
# allows to make a different choice, **in the same format as the file**.
# if more complex behaviour is desired, a separate function is needed.
def GM_adjust_map_core_raw(map_):
    """Makes the necessary changes to the 'raw' input read from core.txt.

    Is expected to not return anything - return value is not caught.

    The core.txt file is stored in Map.rawcore. It has not yet been
    parsed, just loaded into a dictionary. In this dictionary, each
    keyword is its own dictionary key. Most keywords can only occur once
    in the file - those have a list of the 'words' on the line as
    their value. The parameters that are allowed to occur more than once
    have a list as value, in which other lists appear - one for each
    line.

    The purpose of this function is to change this dictionary. Perhaps,
    a rule in core.txt is dependent on a parameter of the map. This
    function can make a decision based on those parameters (stored in
    Map.RunPars).

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    pass


# A function to change the coupling type of an oscillator pair. GMAP can
# only sort oscpairs into the correct couplingmaps based on the name of
# the map of each osc in a pair. However, some maps require different
# coupling maps for different circumstances. This function should return
# the name of the map that should be coupling this pair (instead of itself).
# IF this function does not exist, the map itself is returned by default.
def GM_change_coup_type(map_, system, oscix1, osc1, oscix2, osc2):
    return "DipDip"


# A place to actually do any prepwork. Any preparations should be done here
# (and not pre-frame, for example), as at the time this function is called,
# more information about the oscillator is available (dipole, VEG properties)
def GM_prep_coupling(map_, system, oscixlist, osclist):
    for oscix, osc in zip(oscixlist, osclist):
        # if the map has a function specifically for this map, use it!
        if hasattr(osc.Map.code, "CP_DipDip_calc_dipole"):
            map_.dipole_vec_arr[oscix], dip_pos = (
                osc.Map.code.CP_DipDip_calc_dipole(osc.Map, system, osc))
        else:
            map_.dipole_vec_arr[oscix] = osc.dipole_vec
            dip_pos = osc.dipole_pos

        # save positions in box coordinates
        map_.dipole_pos_arr[oscix] = dip_pos @ system.boxvects_inv


# The GM_calc_coupling function is expected to treat all couplings assigned
# to this map, for a single frame.
def GM_calc_coupling(map_, system, hamiltonian):
    for pair in map_.allpairs:
        oscix1, oscix2 = pair
        J = calc_coupling(
            oscix1, oscix2, map_.dipole_pos_arr, map_.dipole_vec_arr,
            system.boxvects
        )
        hamiltonian[oscix1, oscix2] = J
        hamiltonian[oscix2, oscix1] = J

# Own function of DipDip - njitted for speed.
@njit
def calc_coupling(oscix1, oscix2, pos_arr, vec_arr, boxvects):
    # Used constants:
    # Cm = (1/3.33564) * 10^30 D  (Coulomb meter in Debye)
    # m = 10^10 ang (meter in angstrom)
    # J = 1/hc = (1/1.98644586) * 10^25 1/m
    # => J = 5.03411656 * 10^22 1/cm (joule in wavenumbers)
    # eps_0 = 8.8541878128 F/m = 8.8541878128 C^2/Jm (coulomb squared per
    # joule meter)

    # derived value:
    # 4piEinv = 1/(4 * pi * eps_0) Jm/C^2
    # Gives 5034.11656 cm^-1 * ang*3 Deb^-2

    fourPiEps_inv = np.float32(5034.11656)
    # the positions array is in box-coordinates -> easy subtraction, then
    # move back into cartesian
    d = GM_MF.PBC_back2box(pos_arr[oscix1, :] - pos_arr[oscix2, :], boxvects)
    ir2 = 1/GM_MF.dotprod(d, d)
    ir = np.sqrt(ir2)
    ir3 = ir*ir2
    ir5 = ir3*ir2

    return fourPiEps_inv * (
        GM_MF.dotprod(vec_arr[oscix1], vec_arr[oscix2]) * ir3
        - 3.0 * GM_MF.dotprod(vec_arr[oscix1], d)
        * GM_MF.dotprod(vec_arr[oscix2], d) * ir5)


# A place to do further initialization if a map requires it. Think of
# things like building further lookup tables, for instance.
# (for AmideBB - find neighbours!)
# (Or, for couplings that MUST get information from an oscillator,
# check if that specific function exists)
def GM_post_init(map_, system):
    pass


# A place to do things before the main loop starts (create datastructures
# to be filled in, for example). GEM itself builds the coupling table at
# this point in time. Any preparation stuff that only requires constant
# properties (masses, charges, bonds, for example) should be done here.
def GM_pre_run(map_, system):
    setattr(
        map_, "dipole_vec_arr", np.zeros((system.nosc, 3), dtype="float32"))
    setattr(
        map_, "dipole_pos_arr", np.zeros((system.nosc, 3), dtype="float32"))

    # change dtype of allpair list to suit this map's needs.
    # setattr(map_, "allpairs", np.array(map_.allpairs, dtype='int32').T)
    # setattr(map_, "allpairs_c", np.ctypeslib.as_ctypes(
    #     np.ravel(map_.allpairs)))


# A place to do things before the properties for this frame are being
# calculated. Any preparation stuff that requires frame-dependent
# data should be done here. AIM calculated the CoMs here, GEM also
# builds hamiltonian (as its contents change per frame)
def GM_pre_frame(map_, system):
    pass


# A place to do things with the results from this frame. GEM itself
# writes information like the hamiltonian to files at this point in time.
def GM_post_frame(map_, system):
    pass


# A place to wrap up the entire calculation. GEM itself reports on
# calculation time and treated frames at this point in time.
def GM_post_run(map_, system):
    pass


# A place to define what references to report under what circumstances.
def GM_report_references(map_, system):
    """Returns all references that should be reported for this map.

    This function does not have to account for which outputs are
    actually requested from the program - the text in the reporttext
    field in the references.bib file already does that. It indicates
    for which methods it should be reported, and with which text.

    If this function is absent from a main.py, map_.references will be
    returned in it's entirety. The purpose of this function is to
    return a selection/subset of that dictionary, instead.

    In the case of this mapping, there are different methods for
    computing the different properties of the system, so we only want
    to send those of the selected mapping through, and leave the rest.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)

    Returns
    -------
    references_dict = dict
        This dict should be a slice/part of the full map_.references
        dict. Therefore, an explanation of the map_.references dict:

        The keys in this dictionary are the separate keys listed in the
        reference.bib file's 'mapkey' field. If the reference has
        multiple keys in that field, it will occur multiple times in the
        dictionary, once for each key.

        Associated with each key is a list of
        :class:`~GMAP.src.tools.ReferenceHandler.Reference` objects,
        each of which corresponds to a single entry in the .bib file.
    """

    report_these = {"CP_DipDip": []}

    # loop over all singles maps that have oscillators being coupled by this
    # map
    for singles_map in set(
        [osc.Map for osc in system.oscillators_ordered_coup["DipDip"]]
    ):
        # instead of looking at all references of the map, only look at those
        # that the map itself picked (in case of multiple models and such).
        references = singles_map.code.GM_report_references(singles_map, system)

        # if the map has put aside some references for this map already.
        if "CP_DipDip" in references:
            for reference in references["CP_DipDip"]:
                report_these = add_reference(
                    reference, report_these, singles_map)
        else:
            for mapkey, reflist in references.items():
                for reference in reflist:
                    report_these = add_reference(
                        reference, report_these, singles_map)

    return report_these


# A function specific from this map. Adds some specific references from
# other maps to the dict of references to report for this one.
def add_reference(reference, report_these, singles_map):
    # Now we have a specific reference object. Only if it is meant
    # for dipoles, grab it.
    if "dip" not in reference.reporttext:
        return report_these

    # now, this reference has a dipole reason for being mentioned.
    # So, we'd like to take this reference, but we do make a copy,
    # so we can safely edit the reasons for our own goal.
    new_reference = GM_RH.Reference(reference.input_string)
    new_reference.reporttext = {
        "ham": [
            "dipole moment for the oscillators of type "
            f"{singles_map.name} to be used for the dipole-dipole coupling."
        ]}

    report_these["CP_DipDip"].append(new_reference)
    return report_these
