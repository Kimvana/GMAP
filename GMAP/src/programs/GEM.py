r"""
Usage:

GMAP GEM
GMAP GEM help
    Prints this help.

GMAP GEM demo
    Launches GEM in demo-mode. Performs a basic calculation to
    demonstrate basic use and to verify the program is installed
    correctly.

GMAP GEM run [name of input file] [optional parameters]
    Performs a run of GEM using the parameters specified in the included
    file.


Groningen Electrostatic Maps

The purpose of GEM is to take an MD trajectory and compute the
time-dependent Hamiltonian to be used in electronic spectral calculations.
Instructions on how to deal with specific chromophores have to be included
in the corresponding .emap file.

For more information, check the manual on N/A.
"""


# standard lib imports
import cProfile
import datetime
import subprocess
import sys

# 3rd party lib imports
# import numpy as np

# local imports
import GMAP.src.tools.CLibLoader as GM_CL
import GMAP.src.tools.Exceptions as GM_Ex
import GMAP.src.tools.FileHandler as GM_FH
import GMAP.src.tools.MapReader as GM_MR
import GMAP.src.tools.ParameterParser as GM_PP
import GMAP.src.tools.PhysicsFunctions as GM_PF
import GMAP.src.tools.Plotter as GM_Pl
import GMAP.src.tools.PrintTools as GM_PT
from GMAP.src.tools.PrintTools import devprint as dpr
import GMAP.src.tools.ReferenceHandler as GM_RH
import GMAP.src.tools.SystemReader as GM_SR


# TO DO inside!
def manage_frame(frame, RunPars):
    """Performs all the checks involved with starting a new frame.

    Future/TODO:
    Checks if the new frame should be treated (or is out of range).
    Prints the new frame number, along with an ETA (to know how much
    longer the calculation will take). Also confirms whether there is
    enough time to start on the next batch of frames before time runs
    out.

    Parameters
    ----------
    frame : `MDA.Timestep`
        The frame that will be treated next.
    RunPars : :class:`~GMAP.src.tools.ParameterParser.RunPars`
        The 'main' RunPars instance containing all the basic
        run-defining parameters.
    """

    framenum = frame.frame
    if framenum >= RunPars.stop_frame:
        return True

    # Do a frame number print here! (for ETA type prints)
    relframenum = framenum - RunPars.start_frame

    # If this is the 0th frame, or only first digit is non-zero.
    # We multiply relframenum (the n in nth frame treated) by 10 to make
    # sure that the resulting set is never empty (frames 1 through 9)
    if relframenum == 0 or set(str(relframenum * 10)[1:]) == set("0"):
        # if framenum has form 10^n with n=int
        if str(relframenum)[0] == "1":
            if relframenum != 1:
                GM_PT.Printer.print(2, "")
            verbose_level = 1
        else:
            verbose_level = 2
    else:  # lvl4 prints ETA for each frame.
        verbose_level = 4

    GM_PT.Printer.print(
        4,
        "Current time   | current frame | time elapsed | time to go   | "
        "end time  (est.)"
    )
    print_frame_ETA(
        verbose_level, framenum, RunPars.start_frame, RunPars.stop_frame)

    # Check if there is enough time to do another batch of frames
    return early_stop(framenum, RunPars)


def print_frame_ETA(verbose, framenum, startframe, endframe):
    """Prints some time information about this frame.

    Will report on the time at which the report takes place, the current
    frame, time elapsed, an estimate of the time remaining, and an
    estimate when the program will be done. The following format is
    used:

    Provided format:
    Current time | current frame | time elapsed | time to go   | end time
    Fri 13 HH:MM | xxxyyyzzz     | xxx-xx:xx:xx | xxx-xx:xx:xx | Fri 13 HH:MM

    .. important ::
        These estimates will improve when more frames have already been
        treated. For small systems (proteins, a speed of frames per
        second) any estimate below 10 frames is worthless, after 100
        frames they get usable. For large systems (assemblies, a speed
        of minutes per frame), this estimate will most likely converge
        much faster, but testing is required to know how fast.

        This difference is caused by the contribution of numba jitting.
        This usually takes a few seconds, which is a significant amount
        of time for small systems, but not for larger ones.

    .. note ::
        The weekdays will be reported in the installation(? System?)
        language of the user. As different languages have a shorthand
        for weekdays of a different amount of characters, the program
        has 14 characters reserved (so a few spaces are missing in the
        example above).

    .. note ::
        The estimated time to go (and end time) are based on how long
        earlier frames took. That means that during the first frame
        treated, no estimate can be provided, and wont. The last two
        columns will not be used/filled in on the first frame.

    Parameters
    ----------
    verbose : int
        The verbose level at which the print of this function should be
        performed
    framenum : int
        The frame number at which this function is called
    startframe : int
        The first frame that is treated during the calculation, as
        specified by the user using the parameter start_frame.
    endframe : int
        The (excusive) end point of the calculation, so the first frame
        that won't be treated anymore. As specified by the user using
        the parameter end_frame.
    """

    timer = GM_PT.Printer.Timer
    toprint = []

    # first, add current time (e.g. Fri 13 HH:MM)
    now = datetime.datetime.now()
    datestr = now.strftime("%a %d %H:%M")
    # English has len 12, german has len 11, make it 14 in case any other
    # language needs it... (can't find overview of supported languages)
    toprint.append(f"{datestr: <14}")

    # Then, add current frame number
    toprint.append(f"{framenum: >13}")  # len("currrent frame") == 13

    # Next: time elapsed
    now_ns = timer.get_time("FrameUpdate")
    now_str = GM_PT.time_to_str(now_ns, "s")
    toprint.append(f"{now_str: >12}")  # To fit a max of 999 days.

    if not framenum == startframe:  # if not very first frame of calculation
        # Next: time to go
        start_heavy_ns = timer.get_time("StartLoop")
        ns_per_frame = int((now_ns - start_heavy_ns) / (framenum - startframe))
        ns_to_go = ns_per_frame * (endframe - framenum)
        to_go_str = GM_PT.time_to_str(ns_to_go, "s")
        toprint.append(f"{to_go_str: >12}")  # To fit a max of 999 days.

        # end time
        togo = datetime.timedelta(microseconds=ns_to_go // 1000)
        end_time = now + togo
        datestr = end_time.strftime("%a %d %H:%M")
        # English has len 12, german has len 11, make it 14 in case any other
        # language needs it... (can't find overview of supported languages)
        toprint.append(f"{datestr: <14}")

    GM_PT.Printer.print(verbose, " | ".join(toprint))


def early_stop(framenum, RunPars):
    """Determines whether to stop the calculation early, or to continue.

    This decision is based on the amount of remaining time, used time,
    and frame batch size. Basically, the program divides all frames to
    calculate in batches of a size determined by the user. Every first
    frame of a batch (except the very first batch), the program sees how
    long batches have taken until now, and whether there is enough time
    to finish another.
    If there is not enough time to finish two more, the next batch will
    not start. This is done to ensure that there is also enough time for
    the program to finish things off after the last batch.
    """

    relframenum = framenum - RunPars.start_frame
    if relframenum == 0:  # don't quit on first frame
        return False

    # only consider quitting after completing a batch
    if relframenum % RunPars.batch_size != 0:
        return False

    # now, actually check whether the next batch will fit.
    timer = GM_PT.Printer.Timer
    now_ns = timer.get_time("FrameUpdate")
    start_heavy_ns = timer.get_time("StartLoop")
    ns_per_frame = int((now_ns - start_heavy_ns) / (relframenum))
    avail_time_ns = RunPars.time_limit * 60 * 1000000000

    # if we could do another two batches, allow this batch to continue.
    # why two? because we also need time to finish up the calculation
    # after the last batch.
    if avail_time_ns - now_ns > 2 * ns_per_frame * RunPars.batch_size:
        return False
    else:
        RunPars.end_frame = framenum
        return True


# TO DO inside!
def trj_loop(RunPars, System):
    """Performs the main per-frame loop for GEM.

    Does the last bit of initialization that needs to happen, and then
    treats each frame.

    Parameters
    ----------
    RunPars : :class:`~GMAP.src.tools.ParameterParser.RunPars`
        The 'main' RunPars instance containing all the basic
        run-defining parameters.
    System : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    """

    # first, do precalc
    GM_PT.Printer.add_time(
        3, "Preparing loop over frames", "PrepLoop", "ms")

    # create empty structures, initialize whats needed

    # compare runpar endframe to mda nframes - adjust endframe
    if RunPars.stop_frame >= len(System.universe.trajectory):
        RunPars.stop_frame = len(System.universe.trajectory)

    # Confirm start_frame is still smaller than stop, after the change
    if RunPars.start_frame > RunPars.stop_frame:
        GM_PT.Printer.warning(
            "Encountered an issue with the parameter start_frame. The frame "
            "doesn't exist, as the provided trajectory is too short. In the "
            "current way, there is nothing to do. Please either change the "
            "parameter start_frame to a smaller value, or use a different "
            "trajectory.",
            "SU_WP_17", True, GMAPerrclass=GM_Ex.GmapParameterError
        )

    # let maps prepare for the calculation
    for mapname in System.oscillators_ordered.keys():  # singles
        map_ = RunPars.requested_mapdict[mapname]
        map_.code.GM_pre_run(map_, System)

    # pair maps only need to prepare if couplings are to be calculated.
    if "ham" in RunPars.output_data:
        for mapname in System.oscillators_ordered_coup.keys():  # pairs
            map_ = RunPars.requested_pairmapdict[mapname]
            map_.code.GM_pre_run(map_, System)

    # And in case maps did anything weird...
    RunPars.manage_dtypes()

    # Report on the system we're going to treat.
    GM_FH.write_legend(RunPars, System)

    trj = System.universe.trajectory
    GM_FH.clear_output(RunPars)

    # In case MDA needs a long time to start the loop.
    GM_PT.Printer.add_time(
        3, "Starting loop over frames", "StartLoop", "ms")

    cb = GM_PT.Printer.colors.green_lc
    ct = GM_PT.Printer.colors.clear
    line = f"{cb}════{ct}"
    GM_PT.Printer.print(
        1, f"\n{line} Processing frames {line}", detailed_instructions=[1])
    GM_PT.header(2, "Processing frames", "doublebox_bare")

    # print header for the ETA table (print lvl 4 has header per frame)
    GM_PT.Printer.print(
        1,
        "\nCurrent time   | current frame | time elapsed | time to go   | "
        "end time  (est.)", detailed_instructions=[1, 2, 3]
    )

    for frame in trj[RunPars.start_frame:]:
        GM_PT.Printer.add_time(
            4, "Starting on frame - starting updates", "FrameUpdate", "ms")
        # manage frame number (if not in range, skip, prints, ETA, etc)
        if manage_frame(frame, RunPars):
            break

        # rebuild the frame-specific data (positions, box, etc)
        System.update_properties()
        GM_PT.Printer.add_time(
            4, "done system updates. next: osc updates", "OscUpdate", "ms")
        for oscillator in System.oscillators:
            oscillator.frame_update(System)

        GM_PT.Printer.add_time(
            4, "updates done. next: initialize", "StructInit", "ms")

        # (only if needed) recalc COM

        # initialize output structures (like Ham)
        outputs = GM_PF.generate_output_structures(RunPars, System)

        GM_PT.Printer.add_time(
            4, "initialize done. next: map init", "MapFInit", "ms")

        # call pre-frame funcs of maps
        for mapname in System.oscillators_ordered.keys():  # singles
            map_ = RunPars.requested_mapdict[mapname]
            map_.code.GM_pre_frame(map_, System)

        # pair maps only need to be called if couplings are to be calculated.
        if "ham" in RunPars.output_data:
            for mapname in System.oscillators_ordered_coup.keys():  # pairs
                map_ = RunPars.requested_pairmapdict[mapname]
                map_.code.GM_pre_frame(map_, System)

        GM_PT.Printer.add_time(
            4, "map init done. next: calculation", "Calc", "ms")

        # perform the actual calculations
        outputs = GM_PF.calc_frame(RunPars, System, outputs)

        GM_PT.Printer.add_time(
            4, "calculation done. next: map final", "MapFPost", "ms")

        # call post-frame functions of maps
        for mapname in System.oscillators_ordered.keys():  # singles
            map_ = RunPars.requested_mapdict[mapname]
            map_.code.GM_post_frame(map_, System)

        # pair maps only need to be called if couplings are to be calculated.
        if "ham" in RunPars.output_data:
            for mapname in System.oscillators_ordered_coup.keys():  # pairs
                map_ = RunPars.requested_pairmapdict[mapname]
                map_.code.GM_post_frame(map_, System)

        GM_PT.Printer.add_time(
            4, "map final done. next: write output", "FrameWrite", "ms")

        # write calculated data to files
        GM_FH.write_output(RunPars, frame.frame, outputs)

        GM_PT.Printer.add_time(
            4, "Frame completed. Loading next frame\n", "LoadFrame", "ms")

    cb = GM_PT.Printer.colors.green_lc
    ct = GM_PT.Printer.colors.clear
    GM_PT.Printer.print(
        2, f"\n{cb}====={ct} End of processing frames {cb}====={ct}\n")

    GM_PT.Printer.add_time(
        3, "Frames Completed. Finishing up.", "MapPost", "ms")

    # lastly, do postcalc:
    for mapname in System.oscillators_ordered.keys():  # singles
        map_ = RunPars.requested_mapdict[mapname]
        map_.code.GM_post_run(map_, System)

    # pair maps only need to do postcalc if couplings are to be calculated.
    if "ham" in RunPars.output_data:
        for mapname in System.oscillators_ordered_coup.keys():  # pairs
            map_ = RunPars.requested_pairmapdict[mapname]
            map_.code.GM_post_run(map_, System)

    # print all that the user does not yet know
    # (profiler?)


def print_calculation_summary(RunPars, System):
    """Reports how the calculation went, and some details users might
    want to know.

    Parameters
    ----------
    RunPars : :class:`~GMAP.src.tools.ParameterParser.RunPars`
        The 'main' RunPars instance containing all the basic
        run-defining parameters.
    System : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    """

    # making sure the last 'split' is saved in timer.totals()
    pr = GM_PT.Printer
    GM_PT.Printer.add_time(5, "", "end")

    GM_PT.header(
        1, "Calculation\nsummary", "doublebox_bare", detailed_instructions=[1])
    GM_PT.header(2, "\n  Calculation  \nsummary\n", "doublebox_bare")

    print_time_splits(RunPars)

    print_treated_avail_frames(RunPars, System)

    print_in_output_filenames(RunPars)

    print_relevant_references(RunPars, System)

    end = " ██▓▓▒▒░░"
    start = end[::-1]
    msg = "That was all for today, folks. Thank you, and good night!"
    pr.print(1, f"\n  {start}{msg}{end}")


def print_time_splits(RunPars):
    """Report how much time was spent on what parts of the calculation

    Parameters
    ----------
    RunPars : :class:`~GMAP.src.tools.ParameterParser.RunPars`
        The 'main' RunPars instance containing all the basic
        run-defining parameters.
    """

    def sumavg(*args):
        total_time = pr.Timer.get_total_ns(*args)
        avg_time = total_time // nframes
        tot_str = GM_PT.time_to_str(total_time)
        avg_str = GM_PT.time_to_str(avg_time, "ms")
        return f"{tot_str: >12}  --> {avg_str[-12:]} / frame"

    pr = GM_PT.Printer
    sum_ = pr.Timer.get_total_format
    nframes = RunPars.stop_frame - RunPars.start_frame

    GM_PT.header(2, "Time spent", "doublebox_bare", newlines=(1, 1))
    init_labels = [
        "ParParse", "AddMaps", "MDinit", "MapInit", "ClibLoad", "PrepLoop"]
    f_load = ["StartLoop", "LoadFrame"]
    f_upd = ["FrameUpdate", "PosBox", "COM"]
    f_init = f_upd + ["OscUpdate", "StructInit", "MapFInit"]
    f_calc = ["Calc", "VEGprop", "VEGcalc", "VEGuse"]
    if "ham" in RunPars.output_data:
        f_calc += ["PrepCoup", "CalcCoup"]
    f_post = ["MapFPost", "FrameWrite"]
    perframe = f_init + f_calc + f_post + f_load
    post_labels = ["MapPost"]
    all_labels = init_labels + perframe + post_labels

    pr.print(1, f"Total time:                   {sum_(*all_labels): >12}")
    pr.print(2, f"  Initialization:             {sum_(*init_labels): >12}")
    pr.print(3, f"    Parsing parameters:       {sum_('ParParse'): >12}")
    pr.print(3, f"    Collecting maps:          {sum_('AddMaps'): >12}")
    pr.print(3, f"    Initializing MD system:   {sum_('MDinit'): >12}")
    pr.print(3, f"    Initializing maps:        {sum_('MapInit'): >12}")
    pr.print(3, f"    Loading C libraries:      {sum_('ClibLoad'): >12}")
    pr.print(2, f"  Treating frames:            {sumavg(*perframe)}")
    pr.print(3, f"    Reading frames:           {sumavg(*f_load)}")
    pr.print(3, f"    Per-frame initialization: {sumavg(*f_init)}")
    pr.print(4, f"      Position/box updates:   {sumavg('PosBox')}")
    pr.print(4, f"      Center of Mass:         {sumavg('COM')}")
    pr.print(4, f"      Oscillator updates:     {sumavg('OscUpdate')}")
    pr.print(4, f"      Structure init.:        {sumavg('PosBox')}")
    pr.print(4, f"      Map initialization:     {sumavg('MapFInit')}")
    pr.print(3, f"    Calculation:              {sumavg(*f_calc)}")
    pr.print(4, f"      Calculating estatics:   {sumavg('VEGcalc')}")
    pr.print(4, f"      SingleMap outputs:      {sumavg('VEGuse')}")
    if "ham" in RunPars.output_data:
        pr.print(4, f"      Coupling preparation:   {sumavg('PrepCoup')}")
        pr.print(4, f"      Coupling calculation:   {sumavg('CalcCoup')}")
    pr.print(3, f"    Frame finalization:       {sumavg(*f_post)}")
    pr.print(4, f"      Map finalization:       {sumavg('MapFPost')}")
    pr.print(4, f"      Writing frames:         {sumavg('FrameWrite')}")
    pr.print(2, f"  Calculation finalization:   {sum_(*post_labels): >12}")


def print_treated_avail_frames(RunPars, System):
    """Report what frames from MD are available, which were requested,
    and which actually calculated.

    For now, there is no difference between the requested and calculated
    frames. This will become relevant when the program can stop early
    due to time constraints.

    Parameters
    ----------
    RunPars : :class:`~GMAP.src.tools.ParameterParser.RunPars`
        The 'main' RunPars instance containing all the basic
        run-defining parameters.
    System : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    """

    pr = GM_PT.Printer

    GM_PT.header(2, "MD frames", "doublebox_bare")
    msg = "Frames treated:     " + " " * 12

    # if the calculation was stopped early, the attribute end_frame exists.
    last_frame = getattr(RunPars, "end_frame", RunPars.stop_frame)
    pr.print(1, f"{msg}{RunPars.start_frame}-{last_frame}")
    msg = "Frames requested:   " + " " * 12
    pr.print(2, f"{msg}{RunPars.start_frame}-{RunPars.stop_frame}")
    msg = "Frames available:   " + " " * 12
    pr.print(3, f"{msg}{0}-{len(System.universe.trajectory)}")


def print_in_output_filenames(RunPars):
    """report which files were used during the calculation

    Parameters
    ----------
    RunPars : :class:`~GMAP.src.tools.ParameterParser.RunPars`
        The 'main' RunPars instance containing all the basic
        run-defining parameters.
    """

    def report_files(RunPars, shorthand, printfname, txtverb, binverb):
        if shorthand in RunPars.output_data:
            fname = getattr(RunPars, f"output_{printfname.lower()}_filename")
            if "txt" in RunPars.output_format:
                temp = fname.parent / f"{fname.name}.txt"
                text = f"{printfname} text file:"
                wrapprint(txtverb, f"{text: <28}{temp}", ps)
            if "bin" in RunPars.output_format:
                temp = fname.parent / f"{fname.name}.bin"
                text = f"{printfname} binary file:"
                wrapprint(binverb, f"{text: <28}{temp}", ps)

    def wrapprint(verbose, message, wrap_preline):
        pr.print(verbose, message, wrap_preline=wrap_preline)

    pr = GM_PT.Printer

    GM_PT.header(2, "Files used", "doublebox_bare")
    cb = GM_PT.Printer.colors.green_lc
    ct = GM_PT.Printer.colors.clear
    line = f"{cb}════{ct}"
    GM_PT.Printer.print(
        1, f"\n{line} Files used {line}", detailed_instructions=[1])
    line = f"{cb}========{ct}"
    files = GM_FH.FileLocations
    ps = "  "  # The string to print as pre-wrap
    pr.print(3, f"{line}  Program files and information {line}")
    wrapprint(3, f"Python installation used:   {sys.executable}", ps)
    wrapprint(3, f"GMAP installation used:     {files.script_dir}", ps)
    wrapprint(3, f"Working directory:          {files.cwd}", ps)
    wrapprint(3, f"Program started at:         {files.now_str}", ps)

    pr.print(2, f"\n{line}  Input files {line}")
    wrapprint(1, f"Command issued:             {files.callcommand}", ps)
    wrapprint(2, f"Default parameter file:     {RunPars.defparfilename}", ps)
    wrapprint(2, f"Input parameter file:       {RunPars.inparfilename}", ps)
    wrapprint(1, f"Topology file analyzed:     {RunPars.topology_file}", ps)
    wrapprint(1, f"Trajectory file analyzed:   {RunPars.trajectory_file}", ps)
    mapdirs = ", ".join([str(direc) for direc in RunPars.map_directory])
    wrapprint(2, f"Map directories used:       {mapdirs}", ps)
    wrapprint(2, f"VEG-library file used:      {RunPars.VEG_clib_file}", ps)

    pr.print(2, f"\n{line}  Output files {line}")
    wrapprint(1, f"Logfile generated:          {RunPars.log_filename}", ps)
    fname = RunPars.output_legend_filename
    wrapprint(2, f"Legend file generated:      {fname}", ps)
    if "ham" in RunPars.output_data:
        fname = RunPars.output_couplingvis_filename
        wrapprint(2, f"Coupling visualization:     {fname}", ps)
    report_files(RunPars, "ham", "Hamiltonian", 2, 2)
    report_files(RunPars, "ene", "Energies", 2, 2)
    report_files(RunPars, "dip", "Dipole", 2, 2)
    report_files(RunPars, "ram", "Raman", 2, 2)
    report_files(RunPars, "pos", "Positions", 2, 2)
    report_files(RunPars, "dbp", "Doublepos", 2, 2)
    if RunPars.profiler:
        fname = RunPars.log_profiling_filename
        pr.print(2, f"profiler output:            {fname}")
    if RunPars.profiler_graph:
        fname = RunPars.log_profiling_graph_filename
        pr.print(2, f"profiler visualization:     {fname}")


def print_relevant_references(RunPars, system):
    """report which references should be cited for this calculation

    Parameters
    ----------
    RunPars : :class:`~GMAP.src.tools.ParameterParser.RunPars`
        The 'main' RunPars instance containing all the basic
        run-defining parameters.
    system : :class:`~GMAP.src.tools.SystemReader.System`
        The class containing all the information on the system of the
        MD trajectory.
    """

    GM_PT.header(2, "References to cite", "doublebox_bare")
    cb = GM_PT.Printer.colors.green_lc
    ct = GM_PT.Printer.colors.clear
    line = f"{cb}════{ct}"
    GM_PT.Printer.print(
        1, f"\n{line} References to cite {line}", detailed_instructions=[1])

    all_references = []
    for singles_map in system.oscillators_ordered.keys():
        singles_map = RunPars.requested_mapdict[singles_map]
        all_references.append(
            singles_map.code.GM_report_references(singles_map, system))

    if "ham" in RunPars.output_data:
        for pairs_map in system.oscillators_ordered_coup.keys():
            pairs_map = RunPars.requested_pairmapdict[pairs_map]
            all_references.append(
                pairs_map.code.GM_report_references(pairs_map, system))

    GM_RH.report_references(RunPars, all_references)


# still a placeholder - this function still has to grow. Should in the
# end manage the different run modes, and probably do nothing else?
# This means, a big decision tree: match job, case x: call func_x,
# case y: call func_y, etc. Now, we're basically only doing 1 kind of job.
def GEM(callcommand):
    GM_PT.Printer.add_time(
        3, "Start Parsing GMAP parameters", "ParParse", "ms")
    # step 1 (is GEM in demo mode? to become: What job do we need to do?)
    if callcommand[1] in ("demo"):
        exp_inpfile = False
    else:
        exp_inpfile = True
    # step 2 (very basic cmd line parse)
    job, in_parfile, argslist = GM_PP.parse_commandline(
        callcommand, alljobs, "GMAP GEM", exp_inpfile, True
    )

    # Parameter parsing
    (
        RunPars, singles_mapdict, pairs_mapdict, CmdPars, InPars, DefPars,
        RefPars
    ) = GM_PP.get_parameters(in_parfile, argslist)

    # If requested, profile the run.
    if RunPars.profiler:
        profile = cProfile.Profile()
        profile.enable()

    GM_PT.Printer.add_time(
        3, "Finished GMAP parameters, start adding maps", "AddMaps", "ms")

    # --- end of SU errors ---

    # Map initialization
    GM_MR.manage_maps_singles(RunPars, singles_mapdict)
    GM_MR.manage_maps_pairs(RunPars, pairs_mapdict)
    GM_PT.Printer.add_time(
        2, "Added all maps, start loading C libraries", "ClibLoad", "ms")

    # initialize C library
    GM_CL.VEG_CLib(RunPars)

    GM_PT.Printer.add_time(
        3, "Libraries loaded, start initializing MD system", "MDinit",
        "ms"
    )

    # Looking at MD system - finding oscillators.
    System = GM_SR.System(RunPars)

    GM_PT.Printer.add_time(
        3, "Initialized MD system, start initializing maps", "MapInit", "ms")

    # GEM is now done - let maps initialize as well
    for mapname in System.oscillators_ordered.keys():  # singles
        map_ = RunPars.requested_mapdict[mapname]
        map_.code.GM_post_init(map_, System)

    # Report on what the system looks like (needs singles mapinit)
    System.print_system(RunPars)

    GM_PT.Printer.add_time(
        3, "Initialization complete, start considering pairs", "MDinit", "ms")

    if "ham" in RunPars.output_data:
        # prepare all pair lookup tables.
        System.order_oscillators_pairs(RunPars)

        # let all coupling maps initialize
        for mapname in System.oscillators_ordered_coup.keys():  # pairs
            map_ = RunPars.requested_pairmapdict[mapname]
            map_.code.GM_post_init(map_, System)

        # Save overview of found coupling maps to file.
        GM_Pl.plot_coupling_choices(RunPars, System)

    # Write output parameter file
    GM_FH.write_parameter_file(
        RefPars, RunPars, System, CmdPars, InPars, DefPars)

    # calculate all (requested) frames
    trj_loop(RunPars, System)

    # finalize profiler
    if RunPars.profiler:
        profile.create_stats()
        profile.dump_stats(RunPars.log_profiling_filename)

    if RunPars.profiler_graph:
        strcommand = [
            "gprof2dot", "-f", "pstats",
            RunPars.log_profiling_filename, "-o",
            RunPars.log_profiling_tempfile]
        subprocess.run(strcommand)

        dpr("running dot")
        strcommand = [
            "dot", "-Tpng", "-o", RunPars.log_profiling_graph_filename,
            RunPars.log_profiling_tempfile]
        subprocess.run(strcommand)

        # remove the tempfile again
        RunPars.log_profiling_tempfile.unlink()

    print_calculation_summary(RunPars, System)


# The jobs that GEM can currently execute.
alljobs = [
    "demo",
    "run"
]


def main(callcommand):
    """Fakes behaviour as if called from __main__.

    During normal operation (user types 'GMAP ...' in the command line),
    this function should never be called. This function replicates the
    'normal' behaviour so partial tests are possible.
    """

    if len(callcommand) == 1:
        print(__doc__)
    else:
        GM_FH.FileLocations()
        GEM(callcommand)


if __name__ == "__main__":
    callcommand = sys.argv
    main(callcommand)
