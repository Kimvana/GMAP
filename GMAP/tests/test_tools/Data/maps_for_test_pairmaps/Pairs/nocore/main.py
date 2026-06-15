
def GM_calc_coupling(map_, system, hamiltonian):
    for pair in map_.allpairs:
        oscix1, oscix2 = pair
        J = 1.234
        hamiltonian[oscix1, oscix2] = J
        hamiltonian[oscix2, oscix1] = J
