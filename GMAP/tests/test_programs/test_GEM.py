
# local imports
from GMAP.src.tools import CmdInterface as GM_CI


# Now, we also test whether the program can actually run.
def test_if_runs():
    GM_CI.cmd_interface(["GMAP", "GEM", "run", "../test_inpar.txt"])


def test_early_quit(tmp_path):
    with open(tmp_path / "input_parameters.txt", "w") as fhand:
        fhand.write("output_directory  out\n")
        fhand.write("log_directory   out\n")
    (tmp_path / "out").mkdir()
    GM_CI.cmd_interface([
        "GMAP", "GEM", "run",
        str((tmp_path / "input_parameters.txt").resolve()),
        "--number_frames", "25",
        "--batch_size", "10",
        "--time_limit", "0"
    ])

    with open(tmp_path / "out" / "log.log", encoding="utf-8") as fhand:
        for line in fhand:
            line = line.strip()
            if line.startswith("Frames treated"):
                assert line.split()[2] == "0-10"
            if line.startswith("Frames requested"):
                assert line.split()[2] == "0-25"
            if line.startswith("Frames available"):
                assert line.split()[2] == "0-51"
