r"""
Usage:

GMAP DEPICT
GMAP DEPICT help
    prints this help

GMAP DEPICT calculate [name of input file] [optional parameters]
    Calculates all datapoints for a potential vs estatic_range graph.

GMAP DEPICT show [name of input file] [optional parameters]
    Displays all data calculated previously using calculate.

GMAP DEPICT calcshow [name of input file] [optional parameters]
    Performs the functions of both 'calculate' and 'show'


Dependence of Electrostatic Properties on Individual Charges Taken

The purpose of DEPICT is to visualize how the calculated electrostatic
potential changes with estatic_range. This aids in determining what
value for estatic_range should be used, and what method for calculating
the electrostatic properties.

For more information, check the manual on N/A.
"""


# 3rd party lib imports
import numpy as np
import matplotlib.pyplot as plt

# local imports
import GMAP.src.tools.clib_loader as GM_cl
import GMAP.src.tools.map_reader as GM_mr
import GMAP.src.tools.parameter_parser as GM_pp
import GMAP.src.tools.print_tools as GM_pt
import GMAP.src.tools.system_reader as GM_sr


def calc_data(printer, run_pars, system):
    # GEM is now done - let maps initialize as well
    for mapname in system.oscillators_ordered.keys():
        map_ = run_pars.requested_mapdict[mapname]
        map_.code.GM_pre_run(printer, map_, system)

    # And in case maps did anything weird...
    run_pars.manage_dtypes()
    printer.add_time(3, "Starting on calculation", "ms")

    system.update_properties(run_pars)
    printer.add_time(4, "done system updates. next: osc updates", "ms")

    # only consider a single oscillator
    system.oscillators = [system.oscillators[0]]
    system.nosc = np.int32(1)
    for oscillator in system.oscillators:
        oscillator.frame_update(system)

    printer.add_time(4, "updates done. next: initialize", "ms")

    VEGlib = GM_cl.VEG_CLib()
    # call pre-frame funcs of maps
    for mapname in system.oscillators_ordered.keys():
        map_ = run_pars.requested_mapdict[mapname]
        map_.code.GM_pre_frame(printer, map_, system)

    printer.add_time(4, "map init done. next: calculation", "ms")

    estatics = np.zeros((run_pars.number_frames, 4))
    startsize = run_pars.estatic_range
    for add_r_sphere in range(run_pars.number_frames):
        newsize = startsize + add_r_sphere
        run_pars.estatic_range = np.float32(newsize)
        VEGlib.calcPot_perres_mm(system, run_pars, oscillator)
        estatics[add_r_sphere, 0] = newsize
        estatics[add_r_sphere, 1:] = oscillator.VEGout[:3, 0]
    with open(run_pars.output_estatics_filename, "w") as fhand:
        np.savetxt(fhand, estatics)


def show_data(printer, run_pars):
    with open(run_pars.output_estatics_filename, "r") as fhand:
        data = np.loadtxt(fhand)

    plt.plot(data[:, 0], data[:, 2] - data[:, 1])
    plt.show()
    plt.clf()


def DEPICT(callcommand):
    printer = GM_pt.Printer
    alljobs = [
        "calculate",
        "show",
        "calcshow"
    ]

    job, in_parfile, argslist = GM_pp.parse_commandline(
        callcommand, alljobs, "GMAP DEPICT", True, True
    )

    run_pars, mapdict, _, _, _, _ = GM_pp.get_parameters(
        in_parfile, argslist
    )
    printer.add_time(3, "Parsed GMAP parameters", "ms")

    # end of SU errors

    if job in ("calculate", "calcshow"):
        # do the thing
        GM_mr.manage_maps_singles(run_pars, mapdict)
        printer.add_time(3, "Added all maps", "ms")

        # next - MD system!
        system = GM_sr.System(run_pars)
        printer.add_time(3, "Initialized MD system", "ms")

        # GEM is now done - let maps initialize as well
        for mapname in system.oscillators_ordered.keys():
            map_ = run_pars.requested_mapdict[mapname]
            map_.code.GM_post_init(printer, map_, system)
        printer.add_time(2, "Initialization complete", "ms")

        # initialize C library
        GM_cl.VEG_CLib(run_pars)
        calc_data(printer, run_pars, system)

    if job in ("show", "calcshow"):
        # show the thing
        show_data(printer, run_pars)
