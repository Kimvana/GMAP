
# 3rd party lib imports
# from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.MathFunctions as GM_MF
from GMAP.src.tools.PrintTools import devprint as dpr


def GM_prep_coupling(map_, system, oscixlist, osclist):

    for oscix, osc in zip(oscixlist, osclist):
        COvec = GM_MF.PBC_boxdiff_triclin(
            osc.positions_box[1], osc.positions_box[0], system.boxvects)
        COvec /= GM_MF.vec3_len(COvec)

        CNvec = GM_MF.PBC_boxdiff_triclin(
            osc.positions_box[3], osc.positions_box[0], system.boxvects)
        CNvec = GM_MF.project(COvec, CNvec)
        CNvec /= GM_MF.vec3_len(CNvec)
        z = GM_MF.crossprod(COvec, CNvec)
        z /= GM_MF.vec3_len(z)

        if osc.Map.name == "AmideSC" or osc.resnames[1] != "PRO":
            alpha = map_.alpha_gen
            v = map_.v_gen
        else:
            alpha = map_.alpha_pro
            v = map_.v_pro

        rotmat = np.array([COvec, CNvec, z])
        map_.map_tcc_v[oscix] = np.dot(v, rotmat) * alpha


def GM_calc_coupling(map_, system, hamiltonian):
    for pair in map_.allpairs:
        J = calc_coupling(*pair, map_, system)
        hamiltonian[pair[0], pair[1]] = J
        hamiltonian[pair[1], pair[0]] = J


def calc_coupling(oscix1, oscix2, map_, system):
    osc1 = system.oscillators[oscix1]
    osc2 = system.oscillators[oscix2]

    if osc1.Map.name == "AmideSC" or osc1.resnames[1] != "PRO":
        q1 = map_.q_gen
        dq1 = map_.dq_gen
    else:
        q1 = map_.q_pro
        dq1 = map_.dq_pro

    if osc2.Map.name == "AmideSC" or osc2.resnames[1] != "PRO":
        q2 = map_.q_gen
        dq2 = map_.dq_gen
    else:
        q2 = map_.q_pro
        dq2 = map_.dq_pro

    diff = osc1.positions_box[:, None, :] - osc2.positions_box[None, :, :]
    diff -= np.floor(diff + 0.5)
    diff = diff @ system.boxvects
    r2 = np.sum(diff * diff, axis=2)

    r2[r2 < 0.01] = 1
    ir2 = 1 / r2
    ir = np.sqrt(ir2)
    ir3 = ir * ir2
    ir5 = ir3 * ir2

    v1 = map_.map_tcc_v[oscix1]
    v2 = map_.map_tcc_v[oscix2]
    # dpr(oscix1, oscix2, "\n", v1, "\n", v2)
    # dpr("XXXXXXXXXX", np.sum(v2[None, :, :] * diff, axis=2))

    J = ir * dq1[:, None] * dq2[None, :]
    J -= ir3 * (
        dq1[:, None] * q2[None, :] * np.sum(v2[None, :, :] * diff, axis=2)
        - q1[:, None] * dq2[None, :] * np.sum(v1[:, None, :] * diff, axis=2)
        - np.sum(v1[:, None, :] * v2[None, :, :], axis=2)
        * q1[:, None] * q2[None, :])
    J -= (
        3 * ir5 * q1[:, None] * q2[None, :]
        * np.sum(v2[None, :, :] * diff, axis=2)
        * np.sum(v1[:, None, :] * diff, axis=2))
    J = np.sum(J)
    J *= map_.fourPiEps

    return J


def GM_pre_run(map_, system):
    map_.map_tcc_v = np.zeros((system.nosc, 6, 3), dtype="float32")


def GM_post_init(map_, system):
    # read in all parameters from the constants file
    with open(map_.directory / "constants.txt") as fhand:
        # get_next_line just gets the next not-empty line from the file
        map_.fourPiEps = float(get_next_line(fhand).split()[0])
        map_.alpha_gen = float(get_next_line(fhand).split()[0])
        map_.alpha_pro = float(get_next_line(fhand).split()[0])

        map_.q_gen = np.array(
            [float(item) for item in get_next_line(fhand).split()])
        map_.q_pro = np.array(
            [float(item) for item in get_next_line(fhand).split()])
        map_.dq_gen = np.array(
            [float(item) for item in get_next_line(fhand).split()])
        map_.dq_pro = np.array(
            [float(item) for item in get_next_line(fhand).split()])

        map_.v_gen = np.array([
            [float(item) for item in get_next_line(fhand).split()]
            for _ in range(6)])
        map_.v_pro = np.array([
            [float(item) for item in get_next_line(fhand).split()]
            for _ in range(6)])


def get_next_line(fhand):
    line = ""
    while len(line) == 0:
        line = fhand.readline()
        line = line.split("#")[0].strip()
    return line
