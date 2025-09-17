"""Main code file for dipole-dipole map.

Also is the template/explanation for all that a map can hand directly to
GMAP. Of course, more functions are allowed, and these files can import
other custom python files stored in the same directory (or a directory
therein) as this file.
"""

import numpy as np


# A function to adjust the parameters of the map. For some kinds of
# parameter (especially if theres multiple that are linked), the way
# RunPar is built might not be correct. In this function, the user can
# fix that.
def GM_adjust_RunPars(map_):
    """Makes the necessary changes to Map.RunPar.

    Is expected to not return anything - return value is not caught.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    pass


# # what functions (names) a singles-map must contain for this map to work
# # (name should not include CP_coupmapname part)
# def GM_needs_mapfunc(map_):
#     return ["get_scalar2"]


# # what keywords a singles-map's corefile must contain for this map to work
# # (name should not include the coupmapname part)
# def GM_needs_keyword(map_):
#     return ["scalar1"]


# A place to actually do any prepwork. Any preparations should be done here
# (and not pre-frame, for example), as at the time this function is called,
# more information about the oscillator is available (dipole, VEG properties)
def GM_prep_coupling(map_, system, oscixlist, osclist):
    for oscix, osc in zip(oscixlist, osclist):
        map_.scalar_arr1[oscix] = osc.Map.rawcore["NeedsBoth.scalar1"][0]
        map_.scalar_arr2[oscix] = osc.Map.code.CP_NeedsBoth_get_scalar2()


def GM_calc_coupling(map_, system, oscix1, osc1, oscix2, osc2):
    return (
        map_.scalar_arr1[oscix1] * map_.scalar_arr2[oscix1]
        + map_.scalar_arr1[oscix2] * map_.scalar_arr2[oscix2]
    )


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
    setattr(map_, "scalar_arr1", np.zeros((system.nosc, 1), dtype="float32"))
    setattr(map_, "scalar_arr2", np.zeros((system.nosc, 1), dtype="float32"))


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
