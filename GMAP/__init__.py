r"""
Welcome to the G-MAP package. Currently, this is a work in progress. The
following options are available:

GMAP
GMAP help
    prints this help.

GMAP help [program]
    prints the help for the requested program

GMAP GEM
    This program will create the required files for calculating
    vibrational and electronic spectra.

GMAP Setup
    This program will copy the source code and maps to target folder.
"""


from .src.programs import DEPICT
from .src.programs import GEM
from .src.programs import Setup


alltools = {
    "GEM": GEM,
    "Setup": Setup,
    "DEPICT": DEPICT
}
