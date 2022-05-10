# from .AIM import main as AIMmain
# from .GEM import main as GEMmain
# from .help import main as helpmain

# alltools = {
#     "help": helpmain,
#     "AIM": AIMmain,
#     "GEM": GEMmain
# }

from . import help
from . import AIM
from . import GEM

alltools = {
    "help": help,
    "AIM": AIM,
    "GEM": GEM
}
