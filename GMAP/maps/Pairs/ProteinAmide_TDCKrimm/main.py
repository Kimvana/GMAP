
# 3rd party lib imports
from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.MathFunctions as GM_MF
# from GMAP.src.tools.PrintTools import devprint as dpr


def GM_prep_coupling(map_, system, oscixlist, osclist):
    displace = map_.RunPars.displace
    tanangle = map_.RunPars.tanangle

    for oscix, osc in zip(oscixlist, osclist):
        COvec = GM_MF.PBC_boxdiff_triclin(
            osc.positions_box[1], osc.positions_box[0], system.boxvects)
        map_.dipole_pos_arr[oscix] = (
            osc.positions[0] + displace / GM_MF.vec3_len(COvec) * COvec)
        COvec /= GM_MF.vec3_len(COvec)

        CNvec = GM_MF.PBC_boxdiff_triclin(
            osc.positions_box[3], osc.positions_box[0], system.boxvects)
        CNvec = GM_MF.project(COvec, CNvec)
        CNvec /= GM_MF.vec3_len(CNvec)

        dip_vec = COvec + tanangle * CNvec
        dip_vec /= GM_MF.vec3_len(dip_vec)
        map_.dipole_vec_arr[oscix] = dip_vec

    map_.dipole_pos_arr = map_.dipole_pos_arr @ system.boxvects_inv


def GM_calc_coupling(map_, system, hamiltonian):
    for pair in map_.allpairs:
        J = calc_coupling(
            *pair, map_.dipole_pos_arr, map_.dipole_vec_arr, system.boxvects)
        hamiltonian[pair[0], pair[1]] = J
        hamiltonian[pair[1], pair[0]] = J


@njit
# same as DipDip map, as this is just DipDip map with different dipoles
def calc_coupling(oscix1, oscix2, pos_arr, vec_arr, boxvects):
    fourPiEps_inv = np.float32(580)
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


def GM_pre_run(map_, system):
    map_.dipole_vec_arr = np.zeros((system.nosc, 3), dtype="float32")
    map_.dipole_pos_arr = np.zeros((system.nosc, 3), dtype="float32")
