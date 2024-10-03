"""Main code file for AmideSC map.

Also is the template/explanation for all that a map can hand directly to
GMAP. Of course, more functions are allowed, and these files can import
other custom python files stored in the same directory (or a directory
therein) as this file.
"""

# 3rd party imports
import numpy as np


# GMAP imports
# import GMAP.src.tools.MathFunctions as GM_MF
# import GMAP.src.tools.PhysicsFunctions as GM_PF
from GMAP.src.tools.PrintTools import devprint as dpr

# own module imports
import testimport
import AmideSC_code.calculation_methods as MC_CM
import AmideSC_code.parameter_changer as MC_PC
# from .MapCode import parameter_changer as MC_PC
# from .MapCode import calculation_methods as MC_CM

testimport.importtest()


# A function to adjust the choices made in core.txt. Perhaps, based on
# a detected parameter, a different choice is preferred. This function
# allows to make a different choice, **in the same format as the file**.
# if more complex behaviour is desired, a separate function is needed.
def GM_adjust_map_core_raw(Files, map_):
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
    Files : :class:`~GMAP.src.tools.FileHandler.FileLocations`
        Contains all currently known paths and other file-related
        properties.
        Has to be updated after RunPars is finalized.
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    # still to do (for emaps):
    # assume_length_units, VEG_reference

    # still to do:
    # dipoles, doublepos, xyz?? (or fixed across all maps?)

    dpr("entered the AmideSC map core adjustment")
    MC_PC.adjust_map_core_raw(Files, map_)
    dpr("finished map core adjustment")


# A place to do further initialization if a map requires it. Think of
# things like building further lookup tables, for instance.
# (for AmideBB - find neighbours!)
def GM_post_init(files, map_, system):
    # map_.Core.dipole_gas_phase_array = np.array(map_.Core.dipole_gas_phase)
    # if map_.RunPars.dipole_map_choice == "Jansen":
    #     map_.code.GM_calculate_dipole = (
    #         calculation_methods.calc_dipole_Jansen
    if map_.RunPars.dipole_map_choice == "Torii":
        map_.code.GM_calculate_dipole = MC_CM.calc_dipole_Torii
        map_.Core.dipole_gas_phase = np.float32(map_.Core.dipole_gas_phase)


def GM_calculate_raman(Map, Syst, osc):
    """Returns the raman tensor as a length-6 vector: (xx, xy, xz, yy, yz, zz)

    This method can easily be adapted by other maps for working with raman
    tensors. Make sure that rotation_matrix is an orthogonal 3*3 numpy array
    (so, the vectors making it up are orthonormal).
    Then, the raman_tensor_local can be freely chosen.

    As this is so easily adaptable, this could also be made a standard
    built-in function for GMAP. The only reason this is not the case
    currently, is because raman tensors from maps like this are not
    common yet, so a 'usual' way of determining them has not yet been
    created. Maybe, they won't stay of fixed magnitude in local coordinates
    forever, but depend on sth like VEG or atomic distances in the future.
    """

    # the rotation matrix is available as long as the map specifies
    # estatic_choice to be E or G (which is the case here). It is made
    # available immediately at the beginning of the frame.
    COvec = osc.rotation_matrix[0, :]
    CNvec = osc.rotation_matrix[1, :]
    Zvec = osc.rotation_matrix[2, :]

    theta = 34*np.pi/180
    raman_tensor_local = np.diag([20, 4, 1])

    rotation_matrix = np.zeros((3, 3))
    rotation_matrix[0] = np.cos(theta) * COvec - np.sin(theta) * CNvec
    rotation_matrix[1] = np.sin(theta) * COvec + np.cos(theta) * CNvec
    rotation_matrix[2] = Zvec

    raman_tensor_system = (
        rotation_matrix.T @ raman_tensor_local @ rotation_matrix)

    # old (AIM) version:
    # def tp(vect1):  # tensor product
    #     tensor = np.zeros((6), dtype='float32')
    #     tensor[:3] = vect1[0]*vect1
    #     tensor[3:5] = vect1[1]*vect1[1:]
    #     tensor[5] = vect1[2]*vect1[2]
    #     return tensor
    # Rvec = tp(Rtens[0]) * 20 + tp(Rtens[1]) * 4  + tp(Rtens[2])
    # (here, Rtens is what the current version calls rotation_matrix)

    # now, to numpify this, first, redefine tp.
    # def tp(vect1):
    #     return (vect1[:, None] * vect1[None, :])[np.triu_indices(3)]

    # then, we can do the entire array at once:
    # consts = np.array([20, 4, 1])
    # Rvec = (
    #     Rtens[:, :, None] * Rtens[:, None, :] * consts[:, None, None]
    # ).sum(axis=0)[np.triu_indices(3)]

    # in summation notation (forgetting the triu-indices for flattening):
    # with A_ij == A[i, j]
    # Rvec[i, j] = sum{k=1 -> k=3}(Rtens[k, i] * Rtens[k, j] * consts[k])

    # now, is this equivalent to the new method? Lets derive the summation
    # notation for the new method! (R = rotation matrix, A = local raman tens)
    # Assuming A is diagonal (so only a[i, i] exist)
    # Rvec = R.T @ A @ R
    # Rvec[i, j] = sum{k=1 -> k=3}(R.T[i, k] * (A @ R)[k, j])
    #            = sum{k=1 -> k=3}(R[k, i] * A[k, k] * R[k, j])
    # this is the same as the summation for the AIM version!

    # footnote: what is (A @ R)[k, j]?
    # write it out: (A @ R)[i, j] = sum{k=1 -> k=3}(A[i, k] * R[k, j])
    # but, as only k==i exists for A (rest is 0), this becomes:
    # (A @ R)[i, j] = A[i, i] * R[i, j]

    return raman_tensor_system[np.triu_indices(3)]
