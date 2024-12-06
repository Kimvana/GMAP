
# standard lib imports
import ctypes as ct

# gmap imports
import GMAP.src.tools.CodingTools as GM_CT
import GMAP.src.tools.FileHandler as GM_FH


class TCC_Clib(metaclass=GM_CT.Singleton):
    def __init__(self, map_):
        self.clib = ct.CDLL(str(map_.clibfile))

        self.clib.testme.argypes = []
        self.clib.testme.restype = None

    def testme(self):
        self.clib.testme()


def init_map_for_clib(map_, system):

    map_.clibfile = map_.directory / "src"
    map_.clibfile /= "TCC_clib" + GM_FH.FileLocations.clib_extension
    # map_.clib = TCC_Clib(map_)
    # map_.clib.testme()
