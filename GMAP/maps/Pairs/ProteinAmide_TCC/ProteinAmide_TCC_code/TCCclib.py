
# standard lib imports
import ctypes as ct

# gmap imports
import GMAP.src.tools.CodingTools as GM_CT
import GMAP.src.tools.FileHandler as GM_FH


class TCC_Clib(metaclass=GM_CT.Singleton):
    def __init__(self, map_):
        self.clib = ct.CDLL(str(map_.clibfile))

        self.clib.prep_coupling.argtypes = [
            ct.c_int,  # nosc
            ct.c_int,  # noscats
            ct.POINTER(ct.c_int),  # oscixarr
            ct.POINTER(ct.c_int),  # osc_used_ats,
            ct.POINTER(ct.c_float),  # positions_box
            ct.POINTER(ct.c_float),  # boxvects
            ct.POINTER(ct.c_int),  # dopro
            ct.c_float,  # alpha_gen
            ct.c_float,  # alpha_pro
            ct.POINTER(ct.c_float),  # v_gen
            ct.POINTER(ct.c_float),  # v_pro
            ct.POINTER(ct.c_float)  # tcc_v
        ]
        self.clib.prep_coupling.restype = None

        self.clib.calc_coupling.argtypes = [
            ct.c_int,  # npairs
            ct.POINTER(ct.c_int),  # allpairs
            ct.c_int,  # noscats
            ct.POINTER(ct.c_int),  # oscix_to_ix
            ct.POINTER(ct.c_int),  # osc_used_ats
            ct.POINTER(ct.c_float),  # positions_box
            ct.POINTER(ct.c_float),  # boxvects
            ct.POINTER(ct.c_int),  # dopro
            ct.POINTER(ct.c_float),  # q_gen
            ct.POINTER(ct.c_float),  # q_pro
            ct.POINTER(ct.c_float),  # dq_gen
            ct.POINTER(ct.c_float),  # dq_pro
            ct.c_float,  # fourPiEps
            ct.POINTER(ct.c_float),  # tcc_v
            ct.c_int,  # totosc
            ct.POINTER(ct.c_float)  # hamiltonian
        ]
        self.clib.calc_coupling.restype = None

    def prep_coupling(self, map_, system):
        self.clib.prep_coupling(
            map_.nosc,  # nosc
            map_.noscats,  # noscats
            map_.oscixlist_c,  # oscixarr
            map_.all_used_ats_c,  # osc_used_ats
            system.positions_box_c,  # positions_box
            system.boxvects_c,  # boxvects
            map_.dopro_c,  # dopro
            map_.alpha_gen,  # alpha_gen
            map_.alpha_pro,  # alpha_pro
            map_.v_gen_c,  # v_gen
            map_.v_pro_c,  # v_pro
            map_.map_tcc_v_c  # tcc_v
        )

    def calc_coupling(self, map_, system, hamiltonian_c):
        self.clib.calc_coupling(
            map_.n_allpairs,  # npairs
            map_.allpairs_c,  # allpairs
            map_.noscats,  # noscats
            map_.oscix_to_ix_c,  # oscix_to_ix
            map_.all_used_ats_c,  # osc_used_ats
            system.positions_box_c,  # positions_box
            system.boxvects_c,  # boxvects
            map_.dopro_c,  # dopro
            map_.q_gen_c,  # q_gen
            map_.q_pro_c,  # q_pro
            map_.dq_gen_c,  # dq_gen
            map_.dq_pro_c,  # dq_pro
            map_.fourPiEps,  # fourPiEps
            map_.map_tcc_v_c,  # tcc_v
            system.nosc,  # totosc
            hamiltonian_c  # hamiltonian
        )


def init_map_for_clib(map_, system):

    map_.clibfile = map_.directory / "src"
    map_.clibfile /= "TCC_clib" + GM_FH.FileLocations.clib_extension
    map_.clib = TCC_Clib(map_)
