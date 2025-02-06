import numpy as np
import GMAP.src.tools.MathFunctions as GM_MF

ground_state_charges = np.loadtxt("all_charges_Bchl-c_CHELP.dat",skiprows=2,usecols=1)

def GM_get_rotation_matrix(Map, Syst, osc):
    return np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

def GM_pre_frame(map_,system):
    mg_idx = 26
    for oscillator in system.oscillators_ordered["BChl-c"]:
        resnum = system.resnums[oscillator.used_atoms[mg_idx]]
        # if oscillator.oscix < 10:
        #     print(GM_MF.vec3_len(oscillator.positions[mg_idx] - (system.residues.CoM[resnum] + system.boxdims * np.array([1,0,1]))))

        system.residues.CoM[resnum] = oscillator.positions[mg_idx]
        system.charges[oscillator.used_atoms[0]:oscillator.used_atoms[-1]+1] = ground_state_charges

def GM_pre_run(map_, system):
    osc = system.oscillators_ordered["BChl-c"][0]
    print(osc.used_atoms)
